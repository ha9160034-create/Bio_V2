"""
Unit & Integration Tests for Molecular Dynamics Pipeline
========================================================

Verifies strict adherence to:
- `wiki/gromacs-protein-ligand-protocol.md`
- `wiki/open-force-field-and-ligand-parameterization.md`
- `wiki/simulation-and-benchmark-best-practices.md`
- `wiki/mm-pbsa-and-mm-gbsa.md`
- `AGENTS.md` Guardrails

The tests deliberately do NOT run a full MD simulation. The pipeline is
fail-closed: when GROMACS is absent it must raise an error and must not write a
fabricated trajectory.
"""

import os
import sys
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from md.engine import (
    GromacsBackend,
    detect_gromacs,
    get_gmx_version,
    generate_standalone_scripts,
    resolve_water_model,
)
from md.indexing import build_index_groups, parse_gro_atoms, write_index_from_gro
from md.topology import (
    TopologyBuildError,
    TopologyResult,
    parse_pdbqt_models,
    extract_best_pose_from_docking_run,
    assemble_protein_ligand_complex,
    generate_ligand_position_restraints,
    evaluate_cgenff_penalty,
    clean_receptor_for_md,
    sanitize_resname,
    merge_complex_topology,
    combine_gro_files,
    generate_ligand_itp_openff,
    generate_ligand_itp_cgenff,
    build_master_topology,
    ExtractedPose,
    LigandAtom,
)
from md.simulation import (
    SimulationConfig,
    SimulationError,
    BiophysicalMDSimulator,
    generate_mdp_ions,
    generate_mdp_minimization,
    generate_mdp_nve,
    generate_mdp_nvt,
    generate_mdp_npt,
    generate_mdp_production,
    compute_rhombic_dodecahedron_volume,
    parse_gmx_energy_xvg,
    parse_hbond_num_xvg,
)
from md.analysis import (
    parse_xvg_file,
    kabsch_fit,
    calculate_rmsd_from_trajectory,
    calculate_rmsf_from_trajectory,
    load_mmpbsa_results,
    create_rmsd_plotly,
    create_rmsf_plotly,
    create_hbond_count_plotly,
    create_mmpbsa_breakdown_plotly,
    RMSDRecord,
    RMSFRecord,
)


def _gro_line(resnum: int, resname: str, atomname: str, atomnum: int,
              x: float, y: float, z: float) -> str:
    return (
        f"{resnum:5d}{resname:<5s}{atomname:>5s}{atomnum:5d}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
    )


class TestMDEngineAndScripts(unittest.TestCase):
    """Test GROMACS discovery and standalone script generation."""

    def test_gmx_version_query(self):
        info = get_gmx_version()
        self.assertIn("available", info)
        self.assertIn("version", info)
        self.assertIn("backend", info)

    def test_backend_detection_shape(self):
        backend = detect_gromacs()
        if backend is not None:
            self.assertIn(backend.kind, ("native", "wsl"))
            self.assertTrue(backend.gmx_path)

    def test_wsl_path_translation(self):
        if sys.platform != "win32":
            self.skipTest("WSL path translation is Windows-specific")
        backend = GromacsBackend("wsl", ("wsl", "-e", "/usr/bin/gmx"), "/usr/bin/gmx")
        translated = backend.to_backend_path(Path("C:/Users/DELL/data.gro"))
        self.assertEqual(translated, "/mnt/c/Users/DELL/data.gro")

    def test_water_model_aliases(self):
        self.assertEqual(resolve_water_model("TIP3P"), "tip3p")
        self.assertEqual(resolve_water_model("SPC/E"), "spc")
        self.assertEqual(resolve_water_model("TIP4P-Ew"), "tip4pew")

    def test_standalone_scripts_generation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir)
            scripts = generate_standalone_scripts(
                output_dir=out_dir,
                complex_gro="test_complex.gro",
                topol_top="test_topol.top",
                index_ndx="test_index.ndx",
                box_type="dodecahedron",
                box_distance=1.0,
                production_prefix="test_prod",
                ligand_resname="LIG",
            )
            self.assertTrue(scripts["bash"].is_file())
            self.assertTrue(scripts["batch"].is_file())

            content = scripts["bash"].read_text(encoding="utf-8")
            self.assertIn("editconf -f test_complex.gro -o newbox.gro -bt dodecahedron -d 1.0", content)
            self.assertIn("solvate -cp newbox.gro -cs spc216.gro -p test_topol.top -o solv.gro", content)
            self.assertIn("genion -s ions.tpr -o solv_ions.gro", content)
            # Index must be built from the real system, never hard-coded.
            self.assertIn("make_index.py solv_ions.gro test_index.ndx LIG", content)
            self.assertIn("-n test_index.ndx", content)
            self.assertIn("-fit rot+trans", content)
            self.assertIn("deffnm nve", content)


class TestMDTopology(unittest.TestCase):
    """Pose extraction, receptor cleaning and real topology assembly."""

    def setUp(self):
        self.sample_pdbqt = (
            "MODEL 1\n"
            "REMARK VINA RESULT:    -9.100      0.000      0.000\n"
            "HETATM    1  C1  LIG L   1       1.000   2.000   3.000  1.00  0.00     0.000 C \n"
            "HETATM    2  C2  LIG L   1       2.000   2.000   3.000  1.00  0.00     0.000 C \n"
            "HETATM    3  O1  LIG L   1       3.000   2.000   3.000  1.00  0.00     0.000 O \n"
            "HETATM    4  H1  LIG L   1       1.500   2.800   3.000  1.00  0.00     0.000 H \n"
            "ENDMDL\n"
        )

    def test_pdbqt_pose_parsing_and_complex_preview(self):
        poses = parse_pdbqt_models(self.sample_pdbqt)
        self.assertEqual(len(poses), 1)
        self.assertAlmostEqual(poses[0].score, -9.1)
        self.assertEqual(poses[0].heavy_atom_count, 3)

        with tempfile.TemporaryDirectory() as tmpdir:
            rec = Path(tmpdir) / "receptor.pdb"
            rec.write_text(
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
                "ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C\n"
                "ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 20.00           C\n"
                "ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00 20.00           O\n"
                "END\n",
                encoding="utf-8",
            )
            pdb, gro = assemble_protein_ligand_complex(rec, poses[0], Path(tmpdir))
            self.assertTrue(pdb.is_file())
            self.assertTrue(gro.is_file())
            text = pdb.read_text(encoding="utf-8")
            self.assertIn("LIG", text)
            self.assertEqual(sum(1 for l in gro.read_text().splitlines()[2:-1]), 8)

    def test_ligand_position_restraints(self):
        poses = parse_pdbqt_models(self.sample_pdbqt)
        with tempfile.TemporaryDirectory() as tmpdir:
            out = generate_ligand_position_restraints(poses[0].atoms, Path(tmpdir) / "posre_lig.itp")
            content = out.read_text(encoding="utf-8")
            self.assertIn("[ position_restraints ]", content)
            self.assertIn("1000.0     1000.0     1000.0", content)
            restrained = [
                l for l in content.splitlines()
                if "1000.0" in l and not l.strip().startswith(";")
            ]
            self.assertEqual(len(restrained), 3)

    def test_cgenff_penalty_validation(self):
        self.assertEqual(evaluate_cgenff_penalty(5.2).status, "PASS")
        self.assertEqual(evaluate_cgenff_penalty(25.0).status, "WARNING")
        self.assertEqual(evaluate_cgenff_penalty(75.0).status, "FAIL")

    def test_receptor_cleaner_drops_hetero_and_converts_mse(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            src = Path(tmpdir) / "raw.pdb"
            src.write_text(
                "ATOM      1  N   MSE A   1       0.000   0.000   0.000  1.00 20.00           N\n"
                "ATOM      2  SE  MSE A   1       1.000   0.000   0.000  1.00 20.00          SE\n"
                "ATOM      3  CA AALA A   2       2.000   0.000   0.000  1.00 20.00           C\n"
                "ATOM      4  CA BALA A   2       9.000   0.000   0.000  1.00 20.00           C\n"
                "HETATM    5  O   HOH A 201       5.000   5.000   5.000  1.00 20.00           O\n"
                "END\n",
                encoding="utf-8",
            )
            out, warnings = clean_receptor_for_md(src, Path(tmpdir) / "clean.pdb")
            text = out.read_text(encoding="utf-8")
            self.assertIn("MET", text)
            self.assertIn(" SD ", text)
            self.assertNotIn("MSE", text)
            self.assertNotIn("HOH", text)
            self.assertNotIn(" CA BALA", text)
            self.assertTrue(any("HETATM" in w for w in warnings))

    def test_sanitize_resname(self):
        self.assertEqual(sanitize_resname("JZ4 (2-propylphenol)"), "JZ42P")
        self.assertEqual(sanitize_resname(""), "LIG")

    def test_merge_complex_topology(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            top = base / "topol.top"
            top.write_text(
                '#include "amber99sb-ildn.ff/forcefield.itp"\n'
                '[ moleculetype ]\nProtein_chain_A 3\n'
                '#include "amber99sb-ildn.ff/tip3p.itp"\n'
                '#include "amber99sb-ildn.ff/ions.itp"\n'
                '[ system ]\nTest\n'
                '[ molecules ]\nProtein_chain_A 1\n',
                encoding="utf-8",
            )
            itp = base / "LIG_GMX.itp"
            itp.write_text("[ moleculetype ]\nLIG 3\n", encoding="utf-8")
            out = merge_complex_topology(top, itp, "LIG", base / "merged.top")
            text = out.read_text(encoding="utf-8")
            self.assertIn('#include "LIG_GMX.itp"', text)
            self.assertRegex(text, r"LIG\s+1\s*$")
            # Must be inserted before the protein moleculetype, not after the FF include chain.
            self.assertLess(text.index('#include "LIG_GMX.itp"'), text.index("Protein_chain_A 3"))
            self.assertEqual(text.count("Protein_chain_A 1"), 1)

    def test_combine_gro_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            prot = base / "protein.gro"
            prot.write_text(
                "protein\n"
                "    2\n"
                + _gro_line(1, "ALA", "N", 1, 0.0, 0.0, 0.0) + "\n"
                + _gro_line(1, "ALA", "CA", 2, 0.1, 0.0, 0.0) + "\n"
                "   2.0 2.0 2.0\n",
                encoding="utf-8",
            )
            lig = base / "LIG_GMX.gro"
            lig.write_text(
                "ligand\n"
                "    2\n"
                + _gro_line(1, "LIG", "C1", 1, 0.2, 0.2, 0.2) + "\n"
                + _gro_line(1, "LIG", "H1", 2, 0.3, 0.2, 0.2) + "\n"
                "   2.0 2.0 2.0\n",
                encoding="utf-8",
            )
            out = combine_gro_files(prot, lig, base / "complex.gro", "LIG", expected_ligand_atoms=2)
            atoms, box = parse_gro_atoms(out)
            self.assertEqual(len(atoms), 4)
            self.assertEqual(atoms[-1][1], "LIG")
            self.assertIn("2.0", box)

    def test_legacy_fake_generators_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            atoms = parse_pdbqt_models(self.sample_pdbqt)[0].atoms
            with self.assertRaises(TopologyBuildError):
                generate_ligand_itp_openff(atoms, Path(tmpdir) / "openff.itp")
            with self.assertRaises(TopologyBuildError):
                generate_ligand_itp_cgenff(atoms, Path(tmpdir) / "cgenff.itp")
            with self.assertRaises(TopologyBuildError):
                build_master_topology(Path(tmpdir) / "topol.top", "lig.itp")


class TestMDIndexGeneration(unittest.TestCase):
    """Index groups must be derived from the actual solvated system."""

    def _write_solvated_gro(self, path: Path):
        lines = ["solvated system", "    8"]
        lines.append(_gro_line(1, "ALA", "N", 1, 0.0, 0.0, 0.0))
        lines.append(_gro_line(1, "ALA", "CA", 2, 0.1, 0.0, 0.0))
        lines.append(_gro_line(2, "ALA", "H1", 3, 0.2, 0.0, 0.0))
        lines.append(_gro_line(3, "LIG", "C1", 4, 0.3, 0.0, 0.0))
        lines.append(_gro_line(3, "LIG", "H1", 5, 0.4, 0.0, 0.0))
        lines.append(_gro_line(4, "SOL", "OW", 6, 0.5, 0.0, 0.0))
        lines.append(_gro_line(5, "NA", "NA", 7, 0.6, 0.0, 0.0))
        lines.append(_gro_line(6, "CL", "CL", 8, 0.7, 0.0, 0.0))
        lines.append("   3.0 3.0 3.0")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def test_groups_from_actual_system(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gro = Path(tmpdir) / "solv_ions.gro"
            self._write_solvated_gro(gro)
            groups = build_index_groups(gro, ligand_resname="LIG")
            self.assertEqual(groups["Protein"], [1, 2, 3])
            self.assertEqual(groups["Ligand"], [4, 5])
            self.assertEqual(groups["Ligand_Heavy"], [4])
            self.assertEqual(groups["Water_and_ions"], [6, 7, 8])
            self.assertEqual(groups["Protein_Ligand"], [1, 2, 3, 4, 5])
            self.assertEqual(groups["Backbone"], [1, 2])

    def test_write_index_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gro = Path(tmpdir) / "solv_ions.gro"
            self._write_solvated_gro(gro)
            sizes = write_index_from_gro(gro, Path(tmpdir) / "index.ndx", ligand_resname="LIG")
            text = (Path(tmpdir) / "index.ndx").read_text(encoding="utf-8")
            self.assertIn("[ Protein_Ligand ]", text)
            self.assertIn("[ Water_and_ions ]", text)
            self.assertIn("[ Ligand_Heavy ]", text)
            self.assertEqual(sizes["Ligand"], 2)


class TestMDSimulationMDPAndEnsembles(unittest.TestCase):
    """MDP parameters, unit cell calculations, and ensemble validity."""

    def test_rhombic_dodecahedron_volume_savings(self):
        dodec_vol, cube_vol, savings = compute_rhombic_dodecahedron_volume(1.0, 5.0)
        self.assertGreater(savings, 20.0)
        self.assertLess(savings, 35.0)
        self.assertLess(dodec_vol, cube_vol)

    def test_energy_minimization_mdp(self):
        em_charmm = generate_mdp_minimization(force_field="CHARMM36")
        self.assertIn("vdw-modifier    = force-switch", em_charmm)
        self.assertIn("coulombtype     = PME", em_charmm)
        self.assertIn("emtol           = 1000.0", em_charmm)

        em_amber = generate_mdp_minimization(force_field="Amber14SB")
        self.assertIn("vdw-modifier    = none", em_amber)
        self.assertIn("DispCorr        = EnerPres", em_amber)

    def test_nvt_mdp_guardrails(self):
        nvt = generate_mdp_nvt(temperature_k=300.0, duration_ps=100.0, timestep_ps=0.002)
        self.assertIn("tcoupl                  = v-rescale", nvt)
        self.assertIn("tc-grps                 = Protein_Ligand Water_and_ions", nvt)
        self.assertIn("define                  = -DPOSRES -DPOSRES_LIG", nvt)
        self.assertIn("pcoupl                  = no", nvt)

    def test_npt_mdp_guardrails(self):
        npt = generate_mdp_npt(temperature_k=300.0, pressure_bar=1.0, duration_ps=100.0)
        self.assertIn("pcoupl                  = berendsen", npt)
        self.assertIn("tc-grps                 = Protein_Ligand Water_and_ions", npt)
        self.assertIn("define                  = -DPOSRES -DPOSRES_LIG", npt)

    def test_production_replicate_velocity_generation(self):
        mdp = generate_mdp_production(duration_ns=0.1, gen_vel=True, random_seed=1042)
        self.assertIn("continuation            = no", mdp)
        self.assertIn("gen_vel                 = yes", mdp)
        self.assertIn("gen_seed                = 1042", mdp)
        self.assertIn("pcoupl                  = Parrinello-Rahman", mdp)

    def test_nve_mdp_guardrails(self):
        nve = generate_mdp_nve(duration_ps=10.0, timestep_ps=0.002)
        self.assertIn("tcoupl                  = no", nve)
        self.assertIn("pcoupl                  = no", nve)
        self.assertIn("gen_vel                 = no", nve)
        self.assertIn("continuation            = yes", nve)
        self.assertIn("constraints             = h-bonds", nve)
        self.assertIn("nsteps                  = 5000", nve)

    def test_production_mdp_guardrails_and_rejection_of_berendsen(self):
        prod = generate_mdp_production(duration_ns=1.0, temperature_k=300.0, pressure_bar=1.0)
        self.assertIn("pcoupl                  = Parrinello-Rahman", prod)
        self.assertIn("tcoupl                  = v-rescale", prod)
        self.assertNotIn("-DPOSRES", prod)

        with self.assertRaises(ValueError):
            generate_mdp_production(duration_ns=1.0, barostat="berendsen")

    def test_ions_mdp(self):
        ions = generate_mdp_ions()
        self.assertIn("integrator  = steep", ions)
        self.assertIn("coulombtype = cutoff", ions)


class TestMDAnalysisAndBiophysics(unittest.TestCase):
    """Trajectory parsing, RMSD/RMSF, H-bonds, and real MM-PBSA output parsing."""

    def test_kabsch_fit(self):
        import numpy as np
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        fitted = kabsch_fit(coords + 5.0, coords)
        rmsd = np.sqrt(np.mean(np.sum((fitted - coords) ** 2, axis=1)))
        self.assertAlmostEqual(rmsd, 0.0, places=5)

    def test_load_mmpbsa_results_from_real_file(self):
        sample = (
            "gmx_MMPBSA results\n"
            "Delta G binding:  -42.31 +/- 3.21 kcal/mol\n"
            "VDWAALS   -43.18\n"
            "EEL       -12.34\n"
            "EGB        28.55\n"
            "ESURF      -4.23\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            run = Path(tmpdir)
            (run / "FINAL_RESULTS_MMPBSA.dat").write_text(sample, encoding="utf-8")
            res = load_mmpbsa_results(run)
            self.assertIsNotNone(res)
            self.assertAlmostEqual(res.delta_g_bind, -42.31)
            self.assertAlmostEqual(res.delta_g_bind_std, 3.21)
            self.assertAlmostEqual(res.delta_vdw, -43.18)
            self.assertEqual(res.solvation_model, "MM-GBSA (GB model not specified in output)")
            self.assertFalse(res.internal_dielectric_specified)

    def test_load_mmpbsa_returns_none_without_real_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self.assertIsNone(load_mmpbsa_results(Path(tmpdir)))

    def test_energy_and_hbond_xvg_parsers(self):
        energy = (
            '@ s0 legend "Temperature"\n'
            '@ s1 legend "Density"\n'
            "0.0 300.0 1000.0\n"
            "1.0 301.0 1001.0\n"
        )
        hbond = "@ title \"Hydrogen Bonds\"\n0.0 3\n1.0 4\n2.0 5\n"
        with tempfile.TemporaryDirectory() as tmpdir:
            ep = Path(tmpdir) / "energy.xvg"
            hp = Path(tmpdir) / "hbond.xvg"
            ep.write_text(energy, encoding="utf-8")
            hp.write_text(hbond, encoding="utf-8")
            cols = parse_gmx_energy_xvg(ep)
            self.assertEqual(cols["Temperature"], [300.0, 301.0])
            points = parse_hbond_num_xvg(hp)
            self.assertEqual(points[-1], (2.0, 5.0))

    def test_plotly_chart_generators(self):
        rmsd_records = [RMSDRecord(time_ns=0.1 * i, rmsd_angstrom=1.0 + 0.05 * i) for i in range(10)]
        rmsf_records = [RMSFRecord(resnum=i, resname="ALA", rmsf_angstrom=0.8) for i in range(1, 11)]
        hbond_points = [(0.0, 2.0), (1.0, 3.0), (2.0, 3.0)]

        fig_rmsd = create_rmsd_plotly(rmsd_records)
        fig_rmsf = create_rmsf_plotly(rmsf_records)
        fig_hbond = create_hbond_count_plotly(hbond_points)
        self.assertIsNotNone(fig_rmsd)
        self.assertIsNotNone(fig_rmsf)
        self.assertIsNotNone(fig_hbond)

    def test_mmpbsa_breakdown_plot_from_parsed_file(self):
        from md.analysis import MMPBSAResult
        res = MMPBSAResult(
            delta_g_bind=-42.31,
            delta_g_bind_std=3.21,
            delta_vdw=-43.18,
            delta_elec=-12.34,
            delta_polar_solv=28.55,
            delta_nonpolar_solv=-4.23,
            internal_dielectric=2.0,
            solvation_model="MM-GBSA (GB-OBC2)",
        )
        self.assertIsNotNone(create_mmpbsa_breakdown_plotly(res))


class TestFailClosedSimulation(unittest.TestCase):
    """No GROMACS must mean no simulation, not a synthetic trajectory."""

    def test_simulator_refuses_without_gromacs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cfg = SimulationConfig(run_id="fail_closed")
            with mock.patch("md.simulation.detect_gromacs", return_value=None):
                sim = BiophysicalMDSimulator(config=cfg, output_dir=Path(tmpdir))
                topology = TopologyResult(
                    method="test",
                    force_field="amber99sb-ildn",
                    water_model="tip3p",
                    ligand_resname="LIG",
                    complex_gro=Path(tmpdir) / "complex.gro",
                    complex_pdb=Path(tmpdir) / "complex.pdb",
                    ligand_gro=Path(tmpdir) / "LIG_GMX.gro",
                    ligand_itp=Path(tmpdir) / "LIG_GMX.itp",
                    ligand_posre=Path(tmpdir) / "posre_LIG.itp",
                    topol_top=Path(tmpdir) / "topol.top",
                )
                with self.assertRaises(SimulationError):
                    sim.run_pipeline(topology=topology)
            manifest = json.loads((Path(tmpdir) / "md_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "failed")
            self.assertFalse((Path(tmpdir) / "trajectory_pbc_fit.pdb").exists())


if __name__ == "__main__":
    unittest.main()
