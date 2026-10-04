"""Exercise successful, stale, and interrupted source-generation provenance."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "acctools_source_binding", Path(__file__).resolve().parents[1] / "SKYPE.py"
)
skype = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(skype)
provenance = skype.alignment_provenance


class AlignmentSourceBindingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.binary = self.root / "deps/alignasm/build/alignasm"
        self.binary.parent.mkdir(parents=True)
        self.binary.write_text("alignasm binary A")
        self.mapper = self.root / "minimap2"
        self.mapper.write_text("mapper binary A")
        self.fasta = self.root / "assembly.fa"
        self.fasta.write_text(">query\nACGTACGTACGT\n")
        self.reference = self.root / "reference.fa"
        self.reference.write_text(">chr1\nACGTACGTACGT\n")
        self.gap_tool = self.root / "paf_gap_seq.py"
        self.gap_tool.write_text("# gap extraction fixture\n")
        self.prefix = self.root / "assembly"
        self.primary = self.root / "assembly.paf"
        self.alternate = self.root / "assembly.alt.paf"
        self.selected = self.root / "assembly.aln.paf"
        self.binding = Path(str(self.primary) + ".source_binding.json")
        self.pending = Path(str(self.primary) + ".source_binding.pending.json")
        self.calls = []
        self.fail_primary_once = False
        self.change_source = False
        self.change_index = False
        self.replace_primary_after_alternate = False
        self.empty_primary = False

    def mapper_record(self):
        return dict(provenance.content_signature(self.mapper), version="test-1")

    def run_command(self, command, **kwargs):
        self.calls.append(command)
        if command[0] == str(self.mapper):
            if "-d" in command:
                Path(command[command.index("-d") + 1]).write_text(
                    "index of " + self.reference.read_text())
            else:
                output = Path(command[command.index("-o") + 1])
                primary = output == self.primary
                output.write_text("" if primary and self.empty_primary else
                                  "alignment of " + Path(command[command.index("-o") - 1]).read_text())
                if primary and self.fail_primary_once:
                    self.fail_primary_once = False
                    raise subprocess.CalledProcessError(1, command)
                if primary and self.change_source:
                    self.reference.write_text(self.reference.read_text().replace("A", "T"))
                if primary and self.change_index:
                    Path(command[command.index("-o") - 2]).write_text("corrupt index")
                if not primary and self.replace_primary_after_alternate:
                    self.primary.write_text("replacement primary after original mapper output\n")
        elif command[0] == str(self.binary):
            self.selected.write_text(self.primary.read_text() + self.alternate.read_text())
        elif len(command) > 1 and command[1] == str(self.gap_tool):
            Path(command[-1]).write_text(self.fasta.read_text())
        else:
            self.fail(f"Unexpected subprocess: {command}")

    def invoke(self, force=False):
        with patch.object(skype, "SCRIPT_DIR", str(self.root)), \
                patch.object(provenance, "mapper_identity", side_effect=self.mapper_record), \
                patch.object(skype.subprocess, "run", side_effect=self.run_command):
            return skype.run_alignasm(str(self.prefix), 2, str(self.fasta),
                                      str(self.reference), str(self.binary), force)

    def build_index(self):
        with patch.object(provenance, "mapper_identity", side_effect=self.mapper_record), \
                patch.object(provenance.subprocess, "run", side_effect=self.run_command):
            return provenance.ensure_bound_reference_index(
                self.reference, self.root / "index-cache", 2)

    def test_fresh_generation_has_content_bound_index_and_inputs(self):
        self.invoke()
        record = json.loads(self.binding.read_text())
        self.assertEqual(record["schema"], provenance.SOURCE_SCHEMA)
        self.assertEqual(record["inputs_before"], record["inputs_after"])
        self.assertEqual(record["inputs_after"], provenance.source_inputs(self.fasta, self.reference))
        self.assertEqual(record["reference_index_binding"]["reference"], record["inputs_before"]["reference"])
        self.assertEqual(record["outputs"]["primary_paf"], provenance.content_signature(self.primary))
        self.assertEqual(len(record["generation_commands"]), 3)
        self.assertFalse(self.pending.exists())
        calls = len(self.calls)
        self.invoke()
        self.assertEqual(len(self.calls), calls)

    def test_legacy_primary_with_new_alternate_is_not_attested(self):
        self.primary.write_text("legacy alignment\n")
        self.invoke()
        self.assertFalse(self.binding.exists())
        self.assertFalse(self.pending.exists())
        mappings = [c for c in self.calls if c[0] == str(self.mapper) and "-o" in c]
        self.assertEqual(len(mappings), 1)
        self.assertEqual(mappings[0][-1], str(self.alternate))

    def test_failed_primary_is_retried_and_never_bound(self):
        self.fail_primary_once = True
        with self.assertRaises(subprocess.CalledProcessError):
            self.invoke()
        self.assertTrue(self.primary.exists())
        self.assertTrue(self.pending.exists())
        self.assertFalse(self.binding.exists())
        self.invoke()
        self.assertTrue(self.binding.exists())
        primary_commands = [c for c in self.calls if "-o" in c and c[-1] == str(self.primary)]
        self.assertEqual(len(primary_commands), 2)

    def test_source_change_during_mapping_cannot_be_bound(self):
        self.change_source = True
        with self.assertRaisesRegex(RuntimeError, "Assembly or reference changed"):
            self.invoke()
        self.assertFalse(self.binding.exists())
        self.assertTrue(self.pending.exists())

    def test_index_change_during_mapping_cannot_be_bound(self):
        self.change_index = True
        with self.assertRaisesRegex(RuntimeError, "Reference index binding changed"):
            self.invoke()
        self.assertFalse(self.binding.exists())
        self.assertTrue(self.pending.exists())

    def test_generated_primary_replaced_after_alternate_cannot_be_bound(self):
        self.replace_primary_after_alternate = True
        with self.assertRaisesRegex(RuntimeError, "Generated primary or alternate PAF changed"):
            self.invoke()
        self.assertFalse(self.binding.exists())
        self.assertTrue(self.pending.exists())

    def test_force_failure_invalidates_old_source_binding(self):
        self.invoke()
        self.fail_primary_once = True
        with self.assertRaises(subprocess.CalledProcessError):
            self.invoke(force=True)
        self.assertFalse(self.binding.exists())
        self.assertTrue(self.pending.exists())

    def test_empty_primary_gets_empty_alternate_and_valid_binding(self):
        self.empty_primary = True
        self.invoke()
        record = json.loads(self.binding.read_text())
        self.assertEqual(len(record["generation_commands"]), 1)
        self.assertEqual(self.primary.stat().st_size, 0)
        self.assertEqual(self.alternate.stat().st_size, 0)

    def test_cached_paf_is_not_rebound_to_a_replaced_source(self):
        self.invoke()
        previous = self.binding.read_bytes()
        self.fasta.write_text(self.fasta.read_text().replace("A", "T"))
        self.invoke()
        self.assertEqual(self.binding.read_bytes(), previous)
        self.assertNotEqual(json.loads(previous)["inputs_before"]["fasta"],
                            provenance.content_signature(self.fasta))

    def test_same_size_and_mtime_reference_replacement_changes_index(self):
        first = self.build_index()
        stat = self.reference.stat()
        self.reference.write_text(self.reference.read_text().replace("A", "T"))
        os.utime(self.reference, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        second = self.build_index()
        self.assertNotEqual(first["reference"]["sha256"], second["reference"]["sha256"])
        self.assertNotEqual(first["index"]["path"], second["index"]["path"])

    def test_corrupt_index_is_rebuilt_not_relabelled(self):
        first = self.build_index()
        Path(first["index"]["path"]).write_text("corrupt index")
        second = self.build_index()
        self.assertEqual(first["index"], second["index"])
        self.assertEqual(len([c for c in self.calls if "-d" in c]), 2)

    def test_missing_generation_output_link_cannot_attest_cached_index(self):
        first = self.build_index()
        manifest = Path(first["index"]["path"] + ".source_binding.json")
        record = json.loads(manifest.read_text())
        record["generation_output"]["sha256"] = "0" * 64
        manifest.write_text(json.dumps(record))
        second = self.build_index()
        self.assertEqual(second["generation_output"]["sha256"], second["index"]["sha256"])
        self.assertEqual(len([c for c in self.calls if "-d" in c]), 2)

    def test_mapper_content_replacement_changes_index_even_with_same_version(self):
        first = self.build_index()
        self.mapper.write_text("mapper binary B")
        second = self.build_index()
        self.assertNotEqual(first["index"]["path"], second["index"]["path"])


if __name__ == "__main__":
    unittest.main()
