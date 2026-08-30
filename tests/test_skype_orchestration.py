import contextlib
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


PIPELINE_ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "acctools_skype", PIPELINE_ROOT / "SKYPE.py"
)
skype = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(skype)


def reference_bundle(root):
    return skype.ReferenceBundle(
        name="hs1",
        alignasm_ref=str(root / "ref.fa"),
        depth_ref=str(root / "ref.fa"),
        chr_fai=str(root / "ref.fa.fai"),
        tel_bed=str(root / "telomere.bed"),
        rpt_bed=str(root / "repeat.bed"),
        rcs_bed=str(root / "censat.bed"),
        cyt_bed=str(root / "cytoband.bed"),
        ref_stat=None,
    )


class NativeSkypeOrchestrationTests(unittest.TestCase):
    def run_printed_pipeline(self, root, **overrides):
        values = {
            "CELL_LINE": "sample",
            "PREFIX": str(root / "result"),
            "ctg_paf": str(root / "ctg.paf"),
            "ctg_aln_paf": str(root / "ctg.aln.paf"),
            "utg_paf": str(root / "utg.paf"),
            "utg_aln_paf": str(root / "utg.aln.paf"),
            "depth_loc": str(root / "sample.win.stat.gz"),
            "thread": 2,
            "dep_folder": str(root / "deps"),
            "is_progress": False,
            "skype_force": True,
            "graph_depth": 3,
            "print_args": True,
            "reference_bundle": reference_bundle(root),
        }
        values.update(overrides)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            skype.run_skype(**values)
        return output.getvalue()

    def test_native_commands_skip_removed_stages(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = self.run_printed_pipeline(Path(temporary))
        self.assertIn("01_Preprocess_NClose.py", output)
        self.assertIn("10_Graph_Find_Paths.py", output)
        self.assertNotIn("02_Build_Breakend_Graph_Limited.py", output)
        self.assertIn("23_run_nnls.py", output)
        self.assertIn("31_depth_analysis.py", output)
        self.assertNotIn("24_cluster_weight.py", output)
        self.assertNotIn("30_virtual_sky.py", output)

    def test_vcf_input_uses_the_same_single_nnls_route(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = self.run_printed_pipeline(
                root,
                benchmark_vcf_loc=str(root / "input.vcf"),
            )
        self.assertIn("--vcf_input", output)
        self.assertIn("01_Preprocess_NClose.py", output)
        self.assertIn("10_Graph_Find_Paths.py", output)
        self.assertNotIn("02_Build_Breakend_Graph_Limited.py", output)
        self.assertEqual(output.count("23_run_nnls.py"), 1)
        self.assertNotIn("24_cluster_weight.py", output)
        self.assertNotIn("30_virtual_sky.py", output)

    def test_removed_native_stage_restarts_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            for stage in (24, 30):
                with self.subTest(stage=stage):
                    with self.assertRaisesRegex(ValueError, "skype_start_at"):
                        self.run_printed_pipeline(
                            Path(temporary), skype_start_at=stage
                        )

    def test_legacy_uses_only_combined_stage02(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = self.run_printed_pipeline(
                Path(temporary), legacy=True
            )
        self.assertIn("02_Build_Breakend_Graph_Limited.py", output)
        self.assertNotIn("01_Preprocess_NClose.py", output)
        self.assertNotIn("10_Graph_Find_Paths.py", output)

    def test_restart_numbers_are_mode_specific(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaisesRegex(ValueError, "split 01/10"):
                self.run_printed_pipeline(root, skype_start_at=2)
            with self.assertRaisesRegex(ValueError, "legacy 02"):
                self.run_printed_pipeline(
                    root, legacy=True, skype_start_at=10
                )

            stage10 = self.run_printed_pipeline(root, skype_start_at=10)
            self.assertIn("10_Graph_Find_Paths.py", stage10)
            self.assertNotIn("01_Preprocess_NClose.py", stage10)
            legacy02 = self.run_printed_pipeline(
                root, legacy=True, skype_start_at=2
            )
            self.assertIn("02_Build_Breakend_Graph_Limited.py", legacy02)

    def test_restart_reports_missing_prerequisite_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(
                ValueError, "01_nclose_data.pkl"
            ):
                self.run_printed_pipeline(
                    Path(temporary),
                    skype_start_at=10,
                    print_args=False,
                )

    def test_split_options_are_partitioned_by_stage(self):
        stage01, stage10 = skype.split_stage_options(
            "--check_nclose_count --nclose_count_vaf_threshold 0.2 "
            "--vcf_filter_pass PASS . "
            "--debug-force-nclose chr1:123:+ chr2:456:- "
            "--debug_force_nclose chr3:789:- chr4:1000:+ "
            "--add_indel_graph "
            "--limit_combinations limits.json"
        )
        self.assertEqual(
            stage01,
            [
                "--check-nclose-count",
                "--nclose-count-vaf-threshold",
                "0.2",
                "--vcf-filter-pass",
                "PASS",
                ".",
                "--debug-force-nclose",
                "chr1:123:+",
                "chr2:456:-",
                "--debug-force-nclose",
                "chr3:789:-",
                "chr4:1000:+",
            ],
        )
        self.assertEqual(
            stage10,
            [
                "--add-indel-graph",
                "--limit-combinations",
                "limits.json",
            ],
        )

    def test_split_options_reject_unknown_owned_and_duplicate_flags(self):
        for value, message in (
            ("--unknown", "Unknown"),
            ("--alt replacement.paf", "controlled"),
            ("--verbose --verbose", "Duplicate"),
        ):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, message):
                    skype.split_stage_options(value)

    def test_debug_force_nclose_requires_exactly_two_values(self):
        for value in (
            "--debug-force-nclose chr1:123:+",
            "--debug-force-nclose=chr1:123:+ chr2:456:-",
        ):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "requires 2"):
                    skype.split_stage_options(value)

    def test_add_indel_graph_is_sent_only_to_stage10_with_resources(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = self.run_printed_pipeline(
                Path(temporary), option_02="--add_indel_graph"
            )
        command_lines = output.splitlines()
        stage01 = next(
            line for line in command_lines if "01_Preprocess_NClose.py" in line
        )
        stage10 = next(
            line for line in command_lines if "10_Graph_Find_Paths.py" in line
        )
        self.assertNotIn("--add-indel-graph", stage01)
        self.assertIn("--add-indel-graph", stage10)
        self.assertIn("--main-stat-path", stage10)
        self.assertIn("--censat-bed-path", stage10)

    def test_default_keeps_primary_named_ppc_for_downstream_stages(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = self.run_printed_pipeline(Path(temporary))
        self.assertIn("ctg.aln.paf.ppc.paf", output)
        stage01 = next(
            line
            for line in output.splitlines()
            if "01_Preprocess_NClose.py" in line
        )
        self.assertIn("ctg.aln.paf", stage01)
        self.assertIn("--alt", stage01)
        self.assertIn("utg.aln.paf", stage01)

    def test_native_completion_uses_variant_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = root / "result"
            result.mkdir()
            for filename in (
                "total_cov.png",
                "nclose_report.tsv",
                "SV_call_result.vcf",
                "SKYPE_result.bed",
            ):
                (result / filename).touch()
            output = self.run_printed_pipeline(
                root,
                skype_force=False,
            )
        self.assertEqual(output, "")

    def test_vcf_completion_does_not_require_native_bed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = root / "result"
            result.mkdir()
            for filename in (
                "total_cov.png",
                "nclose_report.tsv",
                "SV_benchmark_result.vcf",
            ):
                (result / filename).touch()
            output = self.run_printed_pipeline(
                root,
                skype_force=False,
                benchmark_vcf_loc=str(root / "input.vcf"),
            )
        self.assertEqual(output, "")

    def test_full_assembly_keeps_its_standalone_entrypoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            depth = root / "sample.win.stat.gz"
            depth.touch()
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                skype.run_full_assembly_skype(
                    cell_line="sample",
                    prefix=str(root / "result"),
                    assembly_path=str(root / "assembly.fa"),
                    assembly_paf=str(root / "assembly.aln.paf"),
                    depth_loc=str(depth),
                    thread=2,
                    dep_folder=str(root / "deps"),
                    is_progress=False,
                    skype_force=True,
                    reference_bundle=reference_bundle(root),
                    print_args=True,
                )
        command = output.getvalue()
        self.assertIn("full_assembly_pipeline.py", command)
        self.assertNotIn("23_run_nnls.py", command)
        self.assertNotIn("31_depth_analysis.py", command)

    def test_full_assembly_rejects_legacy_route(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.object(
                skype,
                "resolve_reference_bundle",
                return_value=reference_bundle(root),
            ), self.assertRaisesRegex(ValueError, "--legacy"):
                skype.analysis(
                    "sample",
                    str(root / "work"),
                    str(root / "contig.fa"),
                    str(root / "unitig.fa"),
                    str(root / "depth.win.stat.gz"),
                    1,
                    str(root / "deps"),
                    False,
                    False,
                    True,
                    skype.run_skype,
                    3,
                    full_assembly=str(root / "assembly.fa"),
                    legacy=True,
                )


if __name__ == "__main__":
    unittest.main()
