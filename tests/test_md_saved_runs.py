"""Tests for the standalone MD finalizer and saved-run loading."""

import json
import tempfile
import unittest
from pathlib import Path

from md.finalize_run import finalize
from md.ui import load_saved_md_result, scan_completed_md_runs


class TestMDFinalizer(unittest.TestCase):
    def _make_run(self, base: Path, completed: bool = True) -> Path:
        run_dir = base / "md_test_run"
        run_dir.mkdir()
        (run_dir / "md.mdp").write_text(
            "; test mdp\nnsteps = 250000\ndt = 0.002\n", encoding="utf-8"
        )
        (run_dir / "md.log").write_text(
            "GROMACS version: 2024.3\n", encoding="utf-8"
        )
        if completed:
            (run_dir / "prod_md.tpr").write_bytes(b"tpr")
            (run_dir / "prod_md.xtc").write_bytes(b"xtc")
            (run_dir / "prod_md.log").write_text(
                "GROMACS version: 2024.3\n"
                "Starting mdrun 'test'\n"
                "Finished mdrun on rank 0\n",
                encoding="utf-8",
            )
        return run_dir

    def test_completed_run_writes_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = self._make_run(Path(directory), completed=True)
            manifest = finalize(run_dir)
            self.assertEqual(manifest["status"], "completed")
            self.assertAlmostEqual(manifest["parameters"]["production_duration_ns"], 0.5)
            self.assertEqual(manifest["engine"]["version"], "2024.3")
            on_disk = json.loads((run_dir / "md_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(on_disk["status"], "completed")

    def test_incomplete_run_is_failed(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = self._make_run(Path(directory), completed=False)
            manifest = finalize(run_dir)
            self.assertEqual(manifest["status"], "failed")
            self.assertIn("error", manifest)

    def test_partial_production_is_failed(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = self._make_run(Path(directory), completed=True)
            (run_dir / "prod_md.log").write_text(
                "Starting mdrun 'test'\nstep 2500\n", encoding="utf-8"
            )
            manifest = finalize(run_dir)
            self.assertEqual(manifest["status"], "failed")
            self.assertIn("Finished mdrun", manifest["error"])

    def test_ligand_resname_is_inferred_from_acpype_gro(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = self._make_run(Path(directory), completed=True)
            (run_dir / "ERLOT_GMX.gro").write_text(
                "test\n    3\n    1ERLOT   C1    1   0.000   0.000   0.000\n",
                encoding="utf-8",
            )
            manifest = finalize(run_dir)
            self.assertEqual(manifest["topology"]["ligand_resname"], "ERLOT")

    def test_replicates_are_counted_from_production_tprs(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = self._make_run(Path(directory), completed=True)
            (run_dir / "rep2_prod_md.tpr").write_bytes(b"tpr")
            (run_dir / "rep2_prod_md.xtc").write_bytes(b"xtc")
            (run_dir / "rep3_prod_md.tpr").write_bytes(b"tpr")
            manifest = finalize(run_dir)
            self.assertEqual(manifest["parameters"]["n_replicates"], 3)

    def test_completed_run_drops_stale_error_and_duration(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = self._make_run(Path(directory), completed=True)
            (run_dir / "md_manifest.json").write_text(
                json.dumps(
                    {
                        "status": "failed",
                        "error": "old failure",
                        "parameters": {"production_duration_ns": 99.0},
                    }
                ),
                encoding="utf-8",
            )
            manifest = finalize(run_dir)
            self.assertEqual(manifest["status"], "completed")
            self.assertNotIn("error", manifest)
            self.assertAlmostEqual(
                manifest["parameters"]["production_duration_ns"], 0.5
            )

    def test_saved_run_is_listed_and_loaded(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir = self._make_run(Path(directory), completed=True)
            finalize(run_dir)
            self.assertEqual(scan_completed_md_runs(Path(directory)), [run_dir])
            loaded = load_saved_md_result(run_dir)
            self.assertIsNotNone(loaded)
            self.assertEqual(loaded["manifest"]["status"], "completed")
            self.assertEqual(loaded["pose_meta"]["ligand_name"], "Ligand")


if __name__ == "__main__":
    unittest.main()
