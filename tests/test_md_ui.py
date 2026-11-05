"""AppTest coverage for the MD results tab.

The test builds a synthetic, clearly-labelled result directory (XVG/PDB files)
and injects it through session state.  It does not run GROMACS; it verifies
that the Streamlit results rendering path does not raise and that the honest
MM-PBSA-unavailable state is shown.
"""

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from streamlit.testing.v1 import AppTest

import md.ui as md_ui
from md.engine import GromacsBackend
from md.simulation import SimulationConfig

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = PROJECT_ROOT / "tmp_gmx_test" / "integration"


class TestMDResultsAppTest(unittest.TestCase):
    def test_md_results_tab_renders_without_exception(self):
        required = ["receptor_clean.pdb", "out_seed_42.pdbqt", "JZ4_3d.sdf", "manifest.json"]
        if not all((ASSET_DIR / name).is_file() for name in required):
            self.skipTest("3HTB/JZ4 integration assets are absent")

        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            docking_run = base / "docking_run"
            docking_run.mkdir()
            for name in required:
                shutil.copy2(ASSET_DIR / name, docking_run / name)

            result_dir = base / "md_run"
            result_dir.mkdir()
            (result_dir / "rmsd_ligand.xvg").write_text(
                "# test\n@ s0 legend \"RMSD\"\n0.0 0.10\n0.5 0.12\n1.0 0.11\n",
                encoding="utf-8",
            )
            (result_dir / "rmsf_protein.xvg").write_text(
                "# test\n1 0.05\n2 0.06\n", encoding="utf-8"
            )
            (result_dir / "hbond_num.xvg").write_text(
                "# test\n0.0 1\n0.5 1\n1.0 0\n", encoding="utf-8"
            )
            (result_dir / "trajectory_pbc_fit.pdb").write_text(
                "MODEL 1\n"
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
                "ENDMDL\n",
                encoding="utf-8",
            )

            manifest = {
                "status": "completed",
                "parameters": {"production_duration_ns": 1.0},
                "engine": {"version": "test"},
                "artifacts": {},
                "replicate_statistics": {},
            }
            pose_meta = {
                "target_pdb": "3HTB",
                "ligand_name": "JZ4",
                "best_vina_score": -8.7,
                "heavy_atoms": 10,
                "total_atoms": 10,
                "ligand_sdf_path": str(docking_run / "JZ4_3d.sdf"),
            }

            app = AppTest.from_file(str(PROJECT_ROOT / "app.py"))
            app.session_state["current_md_result"] = {
                "run_id": "apptest_md",
                "run_dir": result_dir,
                "manifest": manifest,
                "pose_meta": pose_meta,
                "config": SimulationConfig(),
                "ligand_resname": "JZ4",
            }
            fake_backend = GromacsBackend(
                "wsl", ("wsl", "-e", "/usr/bin/gmx"), "/usr/bin/gmx", version="test"
            )
            with patch.object(md_ui, "detect_gromacs", return_value=fake_backend), \
                 patch.object(md_ui, "get_gmx_version", return_value={
                     "available": True, "version": "test", "backend": "wsl",
                     "path": "/usr/bin/gmx", "gpu_acceleration": False,
                 }), \
                 patch.object(md_ui, "check_backend_tools", return_value={
                     "acpype": True, "antechamber": True, "obabel": True,
                     "force_fields": ["amber99sb-ildn"],
                 }), \
                 patch.object(md_ui, "scan_available_docking_runs", return_value=[{
                     "run_id": "fake", "ligand_name": "JZ4", "target_id": "3HTB",
                     "top_score": -8.7, "dir": docking_run,
                 }]):
                app.run(timeout=60)

        self.assertFalse(app.exception, [str(item.value) for item in app.exception])
        self.assertTrue(any("نتائج المسار الفعلي" in str(item.value) for item in app.subheader))
        warnings = " ".join(str(item.value) for item in app.warning)
        self.assertIn("MM-PB(GB)SA", warnings)


if __name__ == "__main__":
    unittest.main()
