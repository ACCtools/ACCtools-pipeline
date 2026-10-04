"""Small coherence/failure controls for future canonical hs1 reference builds."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("hs1_test_pipeline", ROOT / "SKYPE.py")
skype = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(skype)
builder = skype.hs1_reference


class Hs1ReferenceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "canonical.fa"
        self.source.write_text(">chr1 description\nacGTACGTACGT\n>chrY\nACGTCAGTAC\n>chrM\nAA\n")
        self.fai = self.root / "model.fa.fai"
        self.fai.write_text("chr1\t12\t0\t80\t81\nchrY\t10\t0\t80\t81\nchrM\t2\t0\t80\t81\n")
        self.recipe_path = self.root / "recipe.json"
        self.recipe = {
            "schema": builder.RECIPE_SCHEMA, "name": "toy_ypar_noM_v1",
            "coordinate_system": "0-based-half-open",
            "source_fasta_sha256": hashlib.sha256(self.source.read_bytes()).hexdigest(),
            "required_contig_lengths": {"chrY": 10, "chrM": 2},
            "exclude_contigs": ["chrM"], "masks": {"chrY": [[0, 2], [8, 10]]},
        }
        self.save_recipe()
        self.output_root = self.root / "new-references"

    def save_recipe(self):
        self.recipe_path.write_text(json.dumps(self.recipe, sort_keys=True))

    def bind_current_source(self):
        self.recipe["source_fasta_sha256"] = hashlib.sha256(self.source.read_bytes()).hexdigest()
        self.save_recipe()

    def prepare(self):
        return builder.prepare_reference(self.source, self.fai, self.output_root, self.recipe_path)

    def assert_no_published_reference(self):
        self.assertEqual(list(self.output_root.rglob(builder.MANIFEST_NAME)), [])
        self.assertEqual(list(self.output_root.rglob(".building.*")), [])

    def bundle(self, prepared):
        return skype.ReferenceBundle("hs1", prepared["fasta"], str(self.source), str(self.fai),
                                     "tel", "rpt", "rcs", "cyt", None,
                                     prepared["cache_namespace"])

    def test_coordinate_preservation_and_mask_boundaries(self):
        before = self.source.read_bytes()
        prepared = self.prepare()
        self.assertEqual(Path(prepared["fasta"]).read_bytes(),
                         b">chr1\nacGTACGTACGT\n>chrY\nNNGTCAGTNN\n")
        self.assertEqual(builder.read_fai_lengths(prepared["fai"]), {"chr1": 12, "chrY": 10})
        with open(prepared["fasta"], "rb") as handle:
            for row in Path(prepared["fai"]).read_text().splitlines():
                name, length, offset, bases, width = row.split("\t")
                handle.seek(int(offset))
                self.assertEqual(len(handle.read(int(length))), int(length))
                self.assertEqual(int(bases) + 1, int(width))
        manifest = json.loads(Path(prepared["manifest"]).read_text())
        self.assertEqual(manifest["contigs"][1]["masked_base_count"], 4)
        self.assertTrue(manifest["contigs"][2]["excluded"])
        self.assertIsNone(manifest["contigs"][2]["output_sequence_sha256"])
        self.assertEqual(self.source.read_bytes(), before)

    def test_reuse_preserves_artifact_bytes_and_mtimes(self):
        prepared = self.prepare()
        before = {name: (Path(prepared[name]).read_bytes(), Path(prepared[name]).stat().st_mtime_ns)
                  for name in ("fasta", "fai", "manifest")}
        self.assertEqual(self.prepare(), prepared)
        for name, expected in before.items():
            self.assertEqual((Path(prepared[name]).read_bytes(), Path(prepared[name]).stat().st_mtime_ns), expected)

    def test_corrupt_fasta_is_preserved_and_refused_even_with_same_size_mtime(self):
        prepared = self.prepare()
        path = Path(prepared["fasta"])
        stat = path.stat()
        changed = path.read_bytes().replace(b"NNG", b"NNT")
        path.write_bytes(changed)
        os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        with self.assertRaisesRegex(RuntimeError, "invalid and was preserved"):
            self.prepare()
        self.assertEqual(path.read_bytes(), changed)

    def test_corrupt_fai_is_preserved_and_refused(self):
        prepared = self.prepare()
        path = Path(prepared["fai"])
        changed = path.read_text().replace("chrY\t10", "chrY\t9")
        path.write_text(changed)
        with self.assertRaisesRegex(RuntimeError, "invalid and was preserved"):
            self.prepare()
        self.assertEqual(path.read_text(), changed)

    def test_manifest_input_tamper_is_refused(self):
        prepared = self.prepare()
        path = Path(prepared["manifest"])
        manifest = json.loads(path.read_text())
        manifest["inputs"]["source_fasta"]["sha256"] = "0" * 64
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(RuntimeError, "invalid and was preserved"):
            self.prepare()

    def test_changed_source_needs_an_explicit_new_recipe_and_separate_version(self):
        prepared = self.prepare()
        original_output = Path(prepared["fasta"]).read_bytes()
        stat = self.source.stat()
        self.source.write_bytes(self.source.read_bytes().replace(b"acGT", b"tcGT"))
        os.utime(self.source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        with self.assertRaisesRegex(ValueError, "does not match"):
            self.prepare()
        self.bind_current_source()
        updated = self.prepare()
        self.assertNotEqual(updated["fasta"], prepared["fasta"])
        self.assertNotEqual(updated["cache_namespace"], prepared["cache_namespace"])
        self.assertEqual(Path(prepared["fasta"]).read_bytes(), original_output)

    def test_invalid_mask_controls(self):
        cases = ([[0, 2], [8, 11]], [[0, 5], [4, 6]], [[8, 10], [0, 2]],
                 [[-1, 2]], [[2, 2]], [[False, 2]])
        for intervals in cases:
            with self.subTest(intervals=intervals):
                self.recipe["masks"]["chrY"] = intervals
                self.save_recipe()
                with self.assertRaises(ValueError):
                    self.prepare()
                self.assert_no_published_reference()

    def test_invalid_recipe_dictionary_controls(self):
        for mutate in (
            lambda: self.recipe.update(coordinate_system="1-based-closed"),
            lambda: self.recipe.update(exclude_contigs=["missing"]),
            lambda: self.recipe.update(exclude_contigs=["chrM", "chrM"]),
            lambda: self.recipe.update(masks={"chrM": [[0, 1]]}),
            lambda: self.recipe.update(name="../outside"),
        ):
            original = json.loads(json.dumps(self.recipe))
            mutate()
            self.save_recipe()
            with self.assertRaises(ValueError):
                self.prepare()
            self.assert_no_published_reference()
            self.recipe = original

    def test_model_fai_controls(self):
        original = self.fai.read_text()
        cases = (original.replace("chrY\t10", "chrY\t9"),
                 original.replace("chrY\t10\t0\t80\t81\n", ""),
                 original + "chrY\t10\t0\t80\t81\n", "chr1\n")
        for value in cases:
            with self.subTest(value=value):
                self.fai.write_text(value)
                with self.assertRaises(ValueError):
                    self.prepare()
                self.assert_no_published_reference()

    def test_malformed_source_controls(self):
        original = self.source.read_bytes()
        cases = (b"ACGT\n" + original, original + b">chrY\nACGTCAGTAC\n",
                 original.replace(b">chrY\nACGTCAGTAC\n", b""),
                 original.replace(b"ACGTCAGTAC", b"ACGTC?GTAC"),
                 original.replace(b"ACGTCAGTAC", b"ACGTCAGTACA"),
                 original.replace(b"ACGTCAGTAC", b""))
        for value in cases:
            with self.subTest(value=value):
                self.source.write_bytes(value)
                self.bind_current_source()
                with self.assertRaises(ValueError):
                    self.prepare()
                self.assert_no_published_reference()

    def test_input_mutation_during_build_cannot_publish_a_manifest(self):
        original_iter = builder.iter_fasta

        def mutate_source(*args):
            yield from original_iter(*args)
            self.source.write_bytes(self.source.read_bytes().replace(b"acGT", b"tcGT"))

        with patch.object(builder, "iter_fasta", side_effect=mutate_source):
            with self.assertRaisesRegex(RuntimeError, "changed during construction"):
                self.prepare()
        self.assert_no_published_reference()

    def test_reference_bundle_uses_new_artifact_and_cache_namespace(self):
        prepared = self.prepare()
        with patch.object(builder, "prepare_reference", return_value=prepared) as prepare:
            bundle = skype.resolve_reference_bundle(self.root, "hs1", force=True)
        self.assertEqual(bundle.alignasm_ref, prepared["fasta"])
        self.assertEqual(bundle.depth_ref, str(self.root / "chm13v2.0.fa"))
        self.assertEqual(prepare.call_args.args[2], str(self.root / "reference_sets"))
        for function in (skype.alignasm_dir_name, skype.skype_dir_name,
                         skype.full_assembly_skype_dir_name):
            self.assertNotEqual(function(bundle), function("hs1"))
            self.assertTrue(function(bundle).endswith(prepared["cache_namespace"]))
        self.assertNotEqual(skype.full_assembly_paf_cache_path(self.source, bundle),
                            skype.full_assembly_paf_cache_path(self.source, "hs1"))
        self.assertEqual(skype.depth_dir_name(bundle), skype.depth_dir_name("hs1"))
        self.assertEqual(skype.alignasm_dir_name("hg38"), "21_alignasm_hg38")

    def test_explicit_legacy_or_differently_bound_output_directory_is_refused(self):
        bundle = self.bundle(self.prepare())
        output = self.root / "results"
        output.mkdir()
        legacy = output / "legacy.pkl"
        legacy.write_bytes(b"existing result")
        with self.assertRaisesRegex(ValueError, "unbound existing results"):
            skype.prepare_reference_output_dir(output, bundle)
        self.assertEqual(legacy.read_bytes(), b"existing result")
        self.assertEqual(len(list(output.iterdir())), 1)
        fresh = self.root / "fresh-results"
        skype.prepare_reference_output_dir(fresh, bundle)
        skype.prepare_reference_output_dir(fresh, bundle)
        bundle.alignment_cache_namespace += ".changed"
        with self.assertRaisesRegex(ValueError, "different reference"):
            skype.prepare_reference_output_dir(fresh, bundle)

    def test_bound_index_is_separate_and_existing_index_is_preserved(self):
        prepared = self.prepare()
        legacy = self.root / "legacy-short-Y.fa"
        legacy.write_text(">chrY\nNNGTCAGTN\n")
        mapper = self.root / "minimap2"
        mapper.write_text("fake mapper for namespace control")
        provenance = skype.alignment_provenance
        identity = dict(provenance.content_signature(mapper), version="toy-1")
        cache = self.root / "indexes"

        def fake_index(command, **kwargs):
            Path(command[command.index("-d") + 1]).write_bytes(Path(command[-1]).read_bytes())

        with patch.object(provenance.subprocess, "run", side_effect=fake_index):
            old = provenance.ensure_bound_reference_index(legacy, cache, 1, mapper=identity)
            before = Path(old["index"]["path"]).read_bytes()
            new = provenance.ensure_bound_reference_index(prepared["fasta"], cache, 1, mapper=identity)
        self.assertNotEqual(old["index"]["path"], new["index"]["path"])
        self.assertEqual(Path(old["index"]["path"]).read_bytes(), before)
        self.assertEqual(new["reference"]["sha256"], hashlib.sha256(Path(prepared["fasta"]).read_bytes()).hexdigest())


if __name__ == "__main__":
    unittest.main()
