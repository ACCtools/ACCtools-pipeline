"""The native alignment cache must follow the executable and PAF inputs."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "acctools_alignasm_cache", Path(__file__).resolve().parents[1] / "SKYPE.py"
)
skype = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(skype)


class NativeAlignasmCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.prefix = self.root / "sample"
        self.binary = self.root / "alignasm"
        self.binary.write_text("binary version A")
        self.raw = self.root / "sample.paf"
        self.raw.write_text("primary alignment A\n")
        self.alt = self.root / "sample.alt.paf"
        self.alt.write_text("alternative alignment A\n")
        self.output = self.root / "sample.aln.paf"
        self.output.write_text("old alignment of unknown provenance\n")
        self.calls = []

    def run_command(self, command, **kwargs):
        self.assertEqual(command[0], str(self.binary))
        self.calls.append(command)
        self.output.write_text(json.dumps({
            "binary": self.binary.read_text(), "raw": self.raw.read_text(),
            "alt": self.alt.read_text(), "command": command,
        }))

    def invoke(self, thread=1):
        with patch.object(skype.subprocess, "run", side_effect=self.run_command):
            return skype.run_alignasm(str(self.prefix), thread, "unused.fa",
                                      "unused.ref.fa", str(self.binary), False)

    def test_unknown_output_recomputed_then_unchanged_inputs_reused(self):
        self.invoke()
        self.assertEqual(len(self.calls), 1)
        self.invoke()
        self.assertEqual(len(self.calls), 1)

    def test_replaced_binary_invalidates_even_when_size_and_mtime_match(self):
        self.invoke()
        stamp = self.binary.stat()
        self.binary.write_text("binary version B")
        os.utime(self.binary, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
        self.invoke()
        self.assertEqual(len(self.calls), 2)
        self.assertIn("version B", self.output.read_text())

    def test_primary_and_alternative_content_changes_invalidate(self):
        self.invoke()
        for index, path in enumerate((self.raw, self.alt), 2):
            stamp = path.stat()
            path.write_text(path.read_text().replace(" A", " B"))
            os.utime(path, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
            self.invoke()
            self.assertEqual(len(self.calls), index)

    def test_exact_command_change_invalidates(self):
        self.invoke(thread=1)
        self.invoke(thread=2)
        self.assertEqual(len(self.calls), 2)
        self.invoke(thread=2)
        self.assertEqual(len(self.calls), 2)

    def test_modified_output_is_not_reused(self):
        self.invoke()
        self.output.write_text("partial or externally overwritten output\n")
        self.invoke()
        self.assertEqual(len(self.calls), 2)

    def test_failure_cannot_leave_a_valid_cache_marker(self):
        self.invoke()
        self.binary.write_text("binary version B")
        with patch.object(skype.subprocess, "run",
                          side_effect=subprocess.CalledProcessError(1, "alignasm")):
            with self.assertRaises(subprocess.CalledProcessError):
                skype.run_alignasm(str(self.prefix), 1, "unused.fa", "unused.ref.fa",
                                   str(self.binary), False)
        self.assertFalse(Path(str(self.output) + ".alignasm.json").exists())
        self.invoke()
        self.assertEqual(len(self.calls), 2)


if __name__ == "__main__":
    unittest.main()
