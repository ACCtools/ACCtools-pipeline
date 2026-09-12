import contextlib
import importlib.util
import io
import subprocess
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
SKYPE_ROOT = PIPELINE_ROOT.parent / "deps" / "SKYPE"
run_subprocess = subprocess.run


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
        dependency = root / "deps"
        dependency.mkdir(exist_ok=True)
        if not (dependency / "SKYPE").exists():
            (dependency / "SKYPE").symlink_to(SKYPE_ROOT, target_is_directory=True)
        values = {
            "CELL_LINE": "sample",
            "PREFIX": str(root / "result"),
            "ctg_paf": str(root / "ctg.paf"),
            "ctg_aln_paf": str(root / "ctg.aln.paf"),
            "utg_paf": str(root / "utg.paf"),
            "utg_aln_paf": str(root / "utg.aln.paf"),
            "unitig_fasta": str(root / "sample.r.fa"),
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

        def capture_pipeline(command, **kwargs):
            result = run_subprocess(command, capture_output=True, text=True, **kwargs)
            print(result.stdout, end="")
            return result

        output = io.StringIO()
        with contextlib.redirect_stdout(output), patch.object(
            skype.subprocess, "run", side_effect=capture_pipeline
        ) as run:
            skype.run_skype(**values)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(Path(run.call_args.args[0][1]).name, "pipeline.py")
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
            for stage in (2, 24, 30):
                with self.subTest(stage=stage):
                    with self.assertRaises(subprocess.CalledProcessError) as failure:
                        self.run_printed_pipeline(
                            Path(temporary), skype_start_at=stage
                        )
                    self.assertIn("skype_start_at", failure.exception.stderr)


    def test_native_restart_numbers(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaises(subprocess.CalledProcessError) as failure:
                self.run_printed_pipeline(root, skype_start_at=2)
            self.assertIn("native 01/10", failure.exception.stderr)
            stages = {
                1: "01_Preprocess_NClose.py",
                10: "10_Graph_Find_Paths.py",
                11: "11_Ref_Outlier_Contig_Modify.py",
                21: "21_run_depth.py",
                22: "22_save_matrix.py",
                23: "23_run_nnls.py",
                31: "31_depth_analysis.py",
            }
            for start in stages:
                with self.subTest(start=start):
                    output = self.run_printed_pipeline(root, skype_start_at=start)
                    for stage, script in stages.items():
                        if stage >= start:
                            self.assertIn(script, output)
                        else:
                            self.assertNotIn(script, output)

    def test_restart_reports_missing_prerequisite_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(subprocess.CalledProcessError) as failure:
                self.run_printed_pipeline(
                    Path(temporary),
                    skype_start_at=10,
                    print_args=False,
                )
            self.assertIn("01_nclose_data.pkl", failure.exception.stderr)
            self.assertEqual(failure.exception.returncode, 1)




    def test_options_are_sent_to_stage01_with_stage10_resources(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = self.run_printed_pipeline(
                Path(temporary), option_skype="--add_indel_graph"
            )
        command_lines = output.splitlines()
        stage01 = next(
            line for line in command_lines if "01_Preprocess_NClose.py" in line
        )
        stage10 = next(
            line for line in command_lines if "10_Graph_Find_Paths.py" in line
        )
        self.assertIn("--option_skype=--add_indel_graph", stage01)
        self.assertNotIn("--option_skype", stage10)
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
        self.assertIn("--censat-endpoints-dir", stage01)
        preparation = next(line for line in output.splitlines()
                           if "censat_endpoints.py" in line)
        self.assertIn("sample.r.fa", preparation)
        self.assertIn("--reference", preparation)
        self.assertNotIn("--force", preparation)

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



if __name__ == "__main__":
    unittest.main()
