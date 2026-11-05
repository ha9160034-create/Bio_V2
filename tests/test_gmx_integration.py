"""Optional real-GROMACS integration tests.

These tests are skipped by default because they require a working GROMACS
installation (native or WSL) and ACPYPE/AmberTools.  Set
``BIO_RUN_GMX_INTEGRATION=1`` to run them.  They deliberately use the small
3HTB/JZ4 assets already present under ``tmp_gmx_test/`` and never fabricate a
trajectory.
"""

import os
import tempfile
import unittest
from pathlib import Path

from md.engine import check_backend_tools, detect_gromacs
from md.simulation import BiophysicalMDSimulator, SimulationConfig, generate_mdp_minimization
from md.topology import (
    build_protein_topology,
    extract_best_pose_from_docking_run,
    parameterize_ligand_acpype,
    prepare_complex_topology,
)
from docking.ligand_prep import prepare_ligand_from_smiles

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = PROJECT_ROOT / "tmp_gmx_test"


def _integration_enabled() -> bool:
    return os.environ.get("BIO_RUN_GMX_INTEGRATION", "").strip() in {"1", "true", "yes"}


@unittest.skipUnless(_integration_enabled(), "set BIO_RUN_GMX_INTEGRATION=1 to run real GROMACS tests")
class TestRealGromacsPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.backend = detect_gromacs(refresh=True)
        if cls.backend is None:
            raise unittest.SkipTest("GROMACS is not installed natively or in WSL")
        cls.tools = check_backend_tools(cls.backend, refresh=True)

    def test_pdb2gmx_and_grompp(self):
        receptor = ASSET_DIR / "3HTB.pdb"
        if not receptor.is_file():
            self.skipTest("3HTB.pdb integration asset is absent")
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            protein = build_protein_topology(
                self.backend,
                receptor,
                out,
                force_field="amber99sb-ildn",
                water_model="tip3p",
            )
            self.assertGreater(protein["protein_atom_count"], 100)
            pdb2gmx_log = (out / "pdb2gmx.log").read_text(encoding="utf-8")
            self.assertIn("-nointer", pdb2gmx_log)
            self.assertIn("-chainsep", pdb2gmx_log)
            self.assertIn("id_or_ter", pdb2gmx_log)
            (out / "em.mdp").write_text(generate_mdp_minimization(), encoding="utf-8", newline="\n")
            from md.engine import run_gmx_command

            run_gmx_command(
                ["grompp", "-f", "em.mdp", "-c", "protein.gro", "-p", "topol.top",
                 "-o", "em.tpr", "-maxwarn", "5"],
                cwd=out,
                log_name="grompp_em",
                timeout=600,
                backend=self.backend,
            )
            self.assertTrue((out / "em.tpr").is_file())

    def test_tip4pew_topology_swap(self):
        receptor = ASSET_DIR / "3HTB.pdb"
        if not receptor.is_file():
            self.skipTest("3HTB.pdb integration asset is absent")
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            protein = build_protein_topology(
                self.backend,
                receptor,
                out,
                force_field="amber99sb-ildn",
                water_model="tip4pew",
            )
            self.assertTrue(protein["water_validation"]["validated"])
            self.assertEqual(protein["water_validation"]["moleculetype"], "SOL")
            top_text = Path(protein["topol_top"]).read_text(encoding="utf-8")
            self.assertIn("tip4pew.itp", top_text)

    def test_acpype_ligand_has_full_parameter_sections(self):
        if not self.tools.get("acpype"):
            self.skipTest("ACPYPE is not installed in the backend")
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            sdf, _, _ = prepare_ligand_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O", "aspirin", out)
            ligand = parameterize_ligand_acpype(self.backend, sdf, out, ligand_resname="ASP", net_charge=-1)
            text = Path(ligand["itp"]).read_text(encoding="utf-8")
            for section in ("[ atomtypes ]", "[ atoms ]", "[ bonds ]", "[ angles ]", "[ dihedrals ]"):
                self.assertIn(section, text)
            self.assertIn("ASP", text)

    def test_full_topology_for_3htb_jz4(self):
        run = ASSET_DIR / "integration"
        if not (run / "out_seed_42.pdbqt").is_file():
            self.skipTest("3HTB/JZ4 docking asset is absent")
        rec, pose, meta = extract_best_pose_from_docking_run(run)
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            result = prepare_complex_topology(
                rec,
                pose,
                out,
                force_field="amber99sb-ildn",
                water_model="tip3p",
                ligand_sdf=Path(meta["ligand_sdf_path"]) if meta.get("ligand_sdf_path") else None,
                ligand_resname="JZ4",
            )
            self.assertGreater(result.protein_atom_count, 100)
            self.assertEqual(result.ligand_atom_count, 22)
            top_text = Path(result.topol_top).read_text(encoding="utf-8")
            self.assertIn("[ molecules ]", top_text)
            self.assertIn("JZ4", top_text)

    def test_retained_monatomic_zinc_topology(self):
        run = ASSET_DIR / "integration"
        if not (run / "out_seed_42.pdbqt").is_file():
            self.skipTest("3HTB/JZ4 docking asset is absent")
        rec, pose, meta = extract_best_pose_from_docking_run(run)
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            receptor_text = Path(rec).read_text(encoding="utf-8").replace(
                "END\n",
                "HETATM 9000 ZN    ZN A 999       5.000   5.000   5.000  1.00 20.00          ZN\nEND\n",
                1,
            )
            receptor = out / "receptor_zn.pdb"
            receptor.write_text(receptor_text, encoding="utf-8")
            result = prepare_complex_topology(
                receptor,
                pose,
                out,
                force_field="amber99sb-ildn",
                water_model="tip3p",
                ligand_sdf=Path(meta["ligand_sdf_path"]) if meta.get("ligand_sdf_path") else None,
                ligand_resname="JZ4",
                retain_monatomic_ions=True,
            )
            self.assertEqual(len(result.retained_ions), 1)
            self.assertEqual(result.retained_ions[0]["resname"], "ZN")
            self.assertIn("ZN", Path(result.topol_top).read_text(encoding="utf-8"))

    def test_short_full_md_pipeline(self):
        if not self.tools.get("acpype"):
            self.skipTest("ACPYPE is not installed in the backend")
        run = ASSET_DIR / "integration"
        if not (run / "out_seed_42.pdbqt").is_file():
            self.skipTest("3HTB/JZ4 docking asset is absent")
        rec, pose, meta = extract_best_pose_from_docking_run(run)
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            topology = prepare_complex_topology(
                rec,
                pose,
                out,
                force_field="amber99sb-ildn",
                water_model="tip3p",
                ligand_sdf=Path(meta["ligand_sdf_path"]) if meta.get("ligand_sdf_path") else None,
                ligand_resname="JZ4",
            )
            cfg = SimulationConfig(
                run_id="integration_short_md",
                em_max_steps=50,
                nve_duration_ps=1.0,
                nvt_duration_ps=0.005,
                npt_duration_ps=0.005,
                prod_duration_ns=0.0005,
                grompp_maxwarn=5,
                n_replicates=1,
            )
            manifest = BiophysicalMDSimulator(cfg, out).run_pipeline(topology)
            self.assertEqual(manifest["status"], "completed")
            self.assertTrue(Path(manifest["artifacts"]["production_tpr"]).is_file() or (out / manifest["artifacts"]["production_tpr"]).is_file())
            self.assertTrue((out / "prod_md.xtc").is_file())
            self.assertTrue((out / "trajectory_pbc_fit.pdb").is_file())
            self.assertTrue((out / "rmsd_ligand.xvg").is_file())
            self.assertTrue((out / "hbond_num.xvg").is_file())
            self.assertIn("Temperature", manifest["energy_summary"])
            self.assertIn("Density", manifest["energy_summary"])
            self.assertTrue(manifest.get("analysis_reference", "").startswith("em.tpr"))
            self.assertGreater(manifest.get("nve_check", {}).get("frames", 0), 2)
            self.assertEqual(manifest.get("engine_mode"), "real_gromacs_execution")


if __name__ == "__main__":
    unittest.main()
