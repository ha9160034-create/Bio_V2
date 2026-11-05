"""MD regression checks based on wiki/gromacs-protein-ligand-protocol.md,
wiki/force-fields.md and wiki/simulation-and-benchmark-best-practices.md.

Synthetic fixtures below check calculations and UI behavior, not real MD.
"""

import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from streamlit.testing.v1 import AppTest

from md.engine import GromacsBackend, GromacsExecutionError, generate_standalone_scripts, resolve_force_field
from md.simulation import BiophysicalMDSimulator, SimulationConfig, SimulationError


class TestMDReviewRegressions(unittest.TestCase):
    def simulator(self, directory):
        return BiophysicalMDSimulator(
            SimulationConfig(), Path(directory),
            GromacsBackend("native", ("gmx",), "gmx", version="test"),
        )

    def test_nve_slope_uses_actual_ps_even_with_irregular_sampling(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            times = [0, 0.02, 0.06, 0.12, 0.2, 0.4]
            (base / "nve_energy.xvg").write_text(
                '@ s0 legend "Total Energy"\n' +
                "".join(f"{t} {-1000 + t}\n" for t in times), encoding="utf-8",
            )
            (base / "solv_ions.gro").write_text(
                "test\n1\n" +
                f"{1:5d}{'ALA':<5s}{'CA':>5s}{1:5d}{0.:8.3f}{0.:8.3f}{0.:8.3f}\n1 1 1\n",
                encoding="utf-8",
            )
            sim = self.simulator(base)
            sim._gmx = Mock(return_value=(0, "", ""))
            result = sim._extract_nve_drift([])
            self.assertAlmostEqual(result["slope_kj_mol_per_ps"], 1.0)
            self.assertAlmostEqual(result["duration_ps"], 0.4)
            self.assertAlmostEqual(result["drift_kj_mol_per_atom_per_ps"], 1.0)

    def test_nve_duplicate_timestamps_do_not_claim_zero_drift(self):
        with tempfile.TemporaryDirectory() as directory:
            sim = self.simulator(directory)
            sim._gmx = Mock(return_value=(0, "", ""))
            (Path(directory) / "nve_energy.xvg").write_text("0 10\n0 20\n", encoding="utf-8")
            warnings = []
            self.assertEqual(sim._extract_nve_drift(warnings), {})
            self.assertTrue(any("timestamps" in item for item in warnings))

    def test_missing_hbonds_are_unavailable_but_real_zero_is_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base / "rmsd.xvg").write_text("0 0.1\n1 0.2\n", encoding="utf-8")
            sim = self.simulator(base)
            result = sim._replicate_summary(1, "test", "rmsd.xvg", "hbond.xvg")
            self.assertIsNone(result["hbond_mean"])
            (base / "hbond.xvg").write_text("0 0\n1 0\n", encoding="utf-8")
            result = sim._replicate_summary(1, "test", "rmsd.xvg", "hbond.xvg")
            self.assertEqual(result["hbond_mean"], 0.0)
            self.assertGreater(result["hbond_frames"], 0)
            with self.assertRaises(SimulationError):
                sim._replicate_summary(1, "test", "missing.xvg", "hbond.xvg")

    def test_force_field_versions_are_not_substituted(self):
        backend = GromacsBackend("native", ("gmx",), "gmx")
        with patch("md.engine.list_force_fields", return_value=["amber14sb", "amber99sb", "amber99sb-ildn"]):
            self.assertEqual(resolve_force_field(backend, "amber14sb"), "amber14sb")
            self.assertEqual(resolve_force_field(backend, "amber99sb"), "amber99sb")
        with patch("md.engine.list_force_fields", return_value=["amber99sb-ildn"]):
            with self.assertRaises(GromacsExecutionError):
                resolve_force_field(backend, "amber14sb")

    def test_exports_preserve_salt_and_custom_index(self):
        with tempfile.TemporaryDirectory() as directory:
            scripts = generate_standalone_scripts(
                Path(directory), index_ndx="custom.ndx", salt_concentration_molar=0.3,
            )
            for path in scripts.values():
                text = path.read_text(encoding="utf-8")
                self.assertIn("-conc 0.3", text)
                commands = [line for line in text.splitlines() if "trjconv -s" in line]
                self.assertEqual(len(commands), 2)
                self.assertTrue(all("-n custom.ndx" in line for line in commands))

    def test_replicate_rmsd_selects_heavy_atoms_and_backbone(self):
        with tempfile.TemporaryDirectory() as directory:
            sim = self.simulator(directory)
            sim._gmx = Mock(return_value=(0, "", ""))
            sim._mdrun = Mock()
            sim._replicate_summary = Mock(return_value={})
            sim._run_extra_replicate(
                2, SimpleNamespace(topol_top=Path("topol.top")),
                {"npt.mdp": Path("npt.mdp"), "md.mdp": Path("md.mdp")},
            )
            rms_calls = [call for call in sim._gmx.call_args_list if call.args[0][0] == "rms"]
            self.assertEqual(len(rms_calls), 1)
            self.assertEqual(rms_calls[0].kwargs["input_text"], "Backbone\nLigand_Heavy\n")

    def test_results_ui_distinguishes_missing_data_from_measured_zero(self):
        import md.ui as ui

        for has_data in (False, True):
            with self.subTest(has_data=has_data), tempfile.TemporaryDirectory() as directory:
                base = Path(directory)
                receptor = base / "receptor.pdb"
                receptor.write_text("END\n", encoding="utf-8")
                if has_data:
                    (base / "rmsd_ligand.xvg").write_text("0 0\n1 0\n", encoding="utf-8")
                    (base / "hbond_num.xvg").write_text("0 0\n1 0\n", encoding="utf-8")
                meta = {"target_pdb": "TEST", "ligand_name": "LIG", "best_vina_score": -1.0,
                        "heavy_atoms": 1, "total_atoms": 2, "ligand_sdf_path": "test.sdf"}
                app = AppTest.from_string("from md.ui import render_md_tab\nrender_md_tab()")
                app.session_state["current_md_result"] = {
                    "run_id": "synthetic_ui_test", "run_dir": base, "pose_meta": meta,
                    "config": SimulationConfig(), "ligand_resname": "LIG",
                    "manifest": {"status": "completed", "parameters": {"production_duration_ns": 1},
                                 "engine": {"version": "test"}, "replicate_statistics": {
                                     "n_replicates": 2, "hbond_count": {"n": 0, "mean": None},
                                 }},
                }
                with ExitStack() as stack:
                    replacements = {
                        "detect_gromacs": GromacsBackend("native", ("gmx",), "gmx"),
                        "get_gmx_version": {"available": True, "version": "test", "path": "gmx"},
                        "check_backend_tools": {"acpype": True, "force_fields": ["amber99sb-ildn"]},
                        "scan_available_docking_runs": [{"run_id": "test", "ligand_name": "LIG",
                                                         "target_id": "TEST", "top_score": -1., "dir": base}],
                        "extract_best_pose_from_docking_run": (receptor, SimpleNamespace(pdb_block=""), meta),
                        "render_mol_3d": None,
                    }
                    for name, value in replacements.items():
                        stack.enter_context(patch.object(ui, name, return_value=value))
                    app.run(timeout=30)
                self.assertFalse(app.exception, [str(item.value) for item in app.exception])
                metrics = {item.label: item for item in app.metric}
                self.assertEqual(metrics["متوسط RMSD للدواء"].value, "0.00 Å" if has_data else "غير متاح")
                self.assertEqual(metrics["متوسط الروابط الهيدروجينية"].value, "0.0" if has_data else "غير متاح")
                if not has_data:
                    self.assertFalse(metrics["متوسط RMSD للدواء"].delta)


if __name__ == "__main__":
    unittest.main()
