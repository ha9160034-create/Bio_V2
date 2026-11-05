"""Regression tests for the independent audit of Bio_V2-main.

Each test targets a defect that was confirmed in the current code (not merely
reported): sequence-identity denominators, the ``CA``/calcium element trap,
missing pH protonation, top-N versus top-1 validation, the tautological domain
check, binding-mode mismatch before ΔScore, symmetry-blind clustering, stale
clinical presets and the unpinned mutation database.
"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from analysis.alignment import calculate_sequence_identity
from analysis.fetch import validate_mutation_database
from docking.comparison import (
    evaluate_binding_mode_compatibility,
    validate_matched_conditions,
)
from docking.ligand_prep import prepare_ligand_from_smiles
from docking.models import DockingPose, PoseCluster
from docking.presets import HUMAN_CLINICAL_PRESETS, validate_preset
from docking.receptor_prep import _infer_element_from_context
from docking.validation import (
    cluster_poses_across_seeds,
    symmetric_pose_rmsd,
    validate_redocking_top_n,
)
from md.engine import solvent_coordinate_file


def _pdbqt_block(coords):
    lines = []
    for index, (x, y, z) in enumerate(coords, start=1):
        lines.append(
            f"ATOM  {index:5d} C{index:<3d} LIG L   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00    +0.000 C"
        )
    return "\n".join(lines) + "\nENDMDL\n"


def _receptor_pdb(offset=0.0, residues=30):
    return "\n".join(
        f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
        f"{i * 3.8 + offset:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        for i in range(1, residues + 1)
    ) + "\nEND\n"


class TestSequenceIdentity(unittest.TestCase):
    def test_gaps_are_excluded_from_the_denominator(self):
        metrics = calculate_sequence_identity("A-C", "ABC")
        self.assertEqual(metrics["alignment_length"], 3)
        self.assertEqual(metrics["aligned_non_gap"], 2)
        self.assertEqual(metrics["gaps"], 1)
        self.assertEqual(metrics["identity_pct"], 100.0)

    def test_mismatch_and_gap_are_reported_separately(self):
        metrics = calculate_sequence_identity("A-CD", "AEF-")
        self.assertEqual(metrics["identical"], 1)
        self.assertEqual(metrics["aligned_non_gap"], 2)
        self.assertEqual(metrics["gaps"], 2)
        self.assertEqual(metrics["identity_pct"], 50.0)


class TestElementInference(unittest.TestCase):
    def test_protein_ca_is_carbon_not_calcium(self):
        self.assertEqual(_infer_element_from_context("ATOM", "CA", "ALA"), "C")

    def test_calcium_ion_is_calcium(self):
        self.assertEqual(_infer_element_from_context("HETATM", "CA", "CA"), "CA")

    def test_heme_iron_is_iron(self):
        self.assertEqual(_infer_element_from_context("HETATM", "FE", "HEM"), "FE")

    def test_ligand_carbon_named_ca_is_carbon(self):
        self.assertEqual(_infer_element_from_context("HETATM", "CA", "LIG"), "C")

    def test_autodock_types_are_mapped_to_elements(self):
        from md.topology import autodock_type_to_element

        self.assertEqual(autodock_type_to_element("OA", "O1"), "O")
        self.assertEqual(autodock_type_to_element("NA", "N1"), "N")
        self.assertEqual(autodock_type_to_element("A", "C1"), "C")
        self.assertEqual(autodock_type_to_element("HD", "H1"), "H")
        self.assertEqual(autodock_type_to_element("CA", "CA"), "Ca")
        self.assertEqual(autodock_type_to_element("C", "CA"), "C")


class TestLigandProtonation(unittest.TestCase):
    def test_carboxylic_acid_is_deprotonated_at_ph_7_4(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, report = prepare_ligand_from_smiles(
                "CC(=O)OC1=CC=CC=C1C(=O)O", "aspirin", Path(tmpdir), ph=7.4
            )
            self.assertTrue(report["protonation_applied"])
            self.assertEqual(report["formal_charge"], -1)
            self.assertIn("[O-]", report["prepared_smiles"])
            self.assertFalse(report["tautomers_enumerated"])

    def test_aliphatic_amine_is_protonated_at_ph_7_4(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, report = prepare_ligand_from_smiles(
                "NCC1=CC=CC=C1", "benzylamine", Path(tmpdir), ph=7.4
            )
            self.assertTrue(report["protonation_applied"])
            self.assertEqual(report["formal_charge"], 1)
            self.assertIn("[NH3+]", report["prepared_smiles"])

    def test_protonation_failure_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch("docking.ligand_prep._apply_ph_protonation", side_effect=RuntimeError("no pKa model")):
                with self.assertRaises(ValueError):
                    prepare_ligand_from_smiles("CCO", "ethanol", Path(tmpdir), ph=7.4)

    def test_pubchem_cid_is_recorded_with_the_docked_protomer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, report = prepare_ligand_from_smiles(
                "CCO", "ethanol", Path(tmpdir), ph=7.4, pubchem_cid=702
            )
            self.assertEqual(report["pubchem_cid"], 702)
            self.assertTrue(report["pubchem_cid_recorded"])
            self.assertEqual(report["ph"], 7.4)


class TestRedockingValidationGate(unittest.TestCase):
    def _objects(self):
        block = _pdbqt_block([(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)])
        pose = DockingPose(
            rank=1, score=-8.5, seed=42, cluster_id=1,
            rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block=block,
        )
        cluster = PoseCluster(
            cluster_id=1, seed_count=3, supporting_seeds=[42, 101, 2024],
            top_score=-8.5, median_score=-8.5, representative_pose=pose,
        )
        return pose, cluster

    def test_top_n_recovery_does_not_validate_the_protocol(self):
        pose, cluster = self._objects()
        with patch(
            "docking.validation.calculate_mapped_redocking_rmsd",
            side_effect=[(3.0, {}), (1.0, {})],
        ):
            report = validate_redocking_top_n("ref", [pose], [cluster], top_n=10)
        self.assertFalse(report["top1_rmsd_pass"])
        self.assertTrue(report["recovery_rmsd_pass"])
        self.assertFalse(report["protocol_validated"])

    def test_top1_with_seed_support_validates_the_protocol(self):
        pose, cluster = self._objects()
        with patch(
            "docking.validation.calculate_mapped_redocking_rmsd",
            side_effect=[(0.8, {}), (0.8, {})],
        ):
            report = validate_redocking_top_n("ref", [pose], [cluster], top_n=10)
        self.assertTrue(report["top1_rmsd_pass"])
        self.assertTrue(report["top1_cluster_pass"])
        self.assertTrue(report["protocol_validated"])


class TestMatchedComparisonGates(unittest.TestCase):
    def test_domain_mismatch_is_rejected(self):
        kwargs = dict(
            healthy_ligand_hash="abc",
            mutant_ligand_hash="abc",
            healthy_waters_policy="removed",
            mutant_waters_policy="removed",
            healthy_seeds=[42, 101, 2024],
            mutant_seeds=[42, 101, 2024],
            healthy_exhaustiveness=8,
            mutant_exhaustiveness=8,
        )
        valid, _ = validate_matched_conditions("Kinase", "Kinase", **kwargs)
        self.assertTrue(valid)
        valid, reason = validate_matched_conditions("Kinase", "SH2", **kwargs)
        self.assertFalse(valid)
        self.assertIn("Domain mismatch", reason)

    def test_binding_mode_gate_accepts_translated_but_identical_pose(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            healthy_pdb = base / "h.pdb"
            mutant_pdb = base / "m.pdb"
            healthy_pdb.write_text(_receptor_pdb(0.0), encoding="utf-8")
            mutant_pdb.write_text(_receptor_pdb(1.0), encoding="utf-8")

            healthy_pose = _pdbqt_block([(58.0, 0.0, 0.0), (60.0, 0.0, 0.0)])
            mutant_pose = _pdbqt_block([(59.0, 0.0, 0.0), (61.0, 0.0, 0.0)])
            compatible, details = evaluate_binding_mode_compatibility(
                healthy_pose, mutant_pose, healthy_pdb, mutant_pdb, ["A:15"], ["A:15"]
            )
            self.assertTrue(compatible)
            self.assertAlmostEqual(details["ligand_pose_rmsd_angstrom"], 0.0, places=2)

    def test_binding_mode_gate_rejects_different_site(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            healthy_pdb = base / "h.pdb"
            mutant_pdb = base / "m.pdb"
            healthy_pdb.write_text(_receptor_pdb(0.0), encoding="utf-8")
            mutant_pdb.write_text(_receptor_pdb(1.0), encoding="utf-8")
            healthy_pose = _pdbqt_block([(58.0, 0.0, 0.0), (60.0, 0.0, 0.0)])
            mutant_pose = _pdbqt_block([(100.0, 20.0, 0.0), (102.0, 20.0, 0.0)])
            compatible, _ = evaluate_binding_mode_compatibility(
                healthy_pose, mutant_pose, healthy_pdb, mutant_pdb, ["A:15"], ["A:15"]
            )
            self.assertFalse(compatible)


class TestSymmetryAwareClustering(unittest.TestCase):
    def test_flipped_benzene_is_one_cluster(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        mol = Chem.RemoveHs(mol)
        conformer = mol.GetConformer()
        points = np.array([
            [conformer.GetAtomPosition(i).x, conformer.GetAtomPosition(i).y, conformer.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ])
        flipped = points.copy()
        flipped[:, 0] *= -1.0
        flipped[:, 1] *= -1.0

        pose_a = DockingPose(rank=1, score=-8.0, seed=42, cluster_id=0,
                             rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block=_pdbqt_block(points))
        pose_b = DockingPose(rank=2, score=-7.9, seed=101, cluster_id=0,
                             rmsd_lb=0.0, rmsd_ub=0.0, pdbqt_block=_pdbqt_block(flipped))

        raw_rmsd = float(np.sqrt(np.mean(np.sum((points - flipped) ** 2, axis=1))))
        self.assertGreater(raw_rmsd, 2.0)
        sym_rmsd, method = symmetric_pose_rmsd(pose_a.pdbqt_block, pose_b.pdbqt_block)
        self.assertLess(sym_rmsd, 2.0)
        self.assertEqual(method, "rdkit_symmetry_aware_calc_rms")

        clusters = cluster_poses_across_seeds([pose_a, pose_b], rmsd_threshold=2.0)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].seed_count, 2)


class TestClinicalPresets(unittest.TestCase):
    def test_all_presets_pass_schema_validation(self):
        for key, preset in HUMAN_CLINICAL_PRESETS.items():
            self.assertEqual(validate_preset(preset), [], key)

    def test_known_bad_references_were_corrected(self):
        covid = HUMAN_CLINICAL_PRESETS["covid_mpro_paxlovid"]
        self.assertEqual(covid["healthy_id"], "7VH8")
        self.assertEqual(covid["mutant_id"], "8D4N")
        self.assertEqual(covid["healthy_ref_ligand"], "4WI")

        cml = HUMAN_CLINICAL_PRESETS["cml_bcr_abl_imatinib"]
        self.assertEqual(cml["mutant_id"], "3QRJ")

        braf = HUMAN_CLINICAL_PRESETS["melanoma_braf_vemurafenib"]
        self.assertEqual(braf["healthy_id"], "1UWH")
        self.assertEqual(braf["mutant_id"], "3OG7")

        sickle = HUMAN_CLINICAL_PRESETS["sickle_cell_voxelotor"]
        self.assertEqual(sickle["healthy_chain"], "A")
        self.assertFalse(sickle["matched_comparison_supported"])

    def test_no_ligand_is_claimed_on_the_apo_t790m_structure(self):
        preset = HUMAN_CLINICAL_PRESETS["egfr_erlotinib_t790m"]
        self.assertEqual(preset["mutant_id"], "2JIT")
        self.assertEqual(preset["mutant_ref_ligand"], "")


class TestMutationDatabaseAndSolvent(unittest.TestCase):
    def test_mutation_schema_rejects_malformed_entries(self):
        self.assertTrue(validate_mutation_database({"1M17": "2JIT"})[0])
        self.assertFalse(validate_mutation_database({"BAD": "2JIT"})[0])
        self.assertFalse(validate_mutation_database(["1M17"])[0])
        self.assertFalse(validate_mutation_database({"1M17": "not-a-pdb"})[0])

    def test_local_mutation_database_is_valid(self):
        path = Path(__file__).resolve().parents[1] / "data" / "mutation.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        valid, message = validate_mutation_database(data)
        self.assertTrue(valid, message)

    def test_solvent_template_matches_water_model(self):
        self.assertEqual(solvent_coordinate_file("tip3p"), "spc216.gro")
        self.assertEqual(solvent_coordinate_file("SPC/E"), "spc216.gro")
        self.assertEqual(solvent_coordinate_file("TIP4P-Ew"), "tip4p.gro")


class TestFailClosedPipeline(unittest.TestCase):
    def test_index_failure_writes_failed_manifest(self):
        from unittest.mock import MagicMock

        from md.engine import GromacsBackend
        from md.simulation import BiophysicalMDSimulator, SimulationConfig, SimulationError
        from md.topology import TopologyResult

        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            for name in ("complex.gro", "complex.pdb", "LIG_GMX.gro",
                         "LIG_GMX.itp", "posre_LIG.itp"):
                (out / name).write_text("placeholder\n", encoding="utf-8")
            (out / "topol.top").write_text(
                '#include "amber99sb-ildn.ff/forcefield.itp"\n'
                "[ molecules ]\nProtein 1\nLIG 1\n",
                encoding="utf-8",
            )
            topology = TopologyResult(
                method="test",
                force_field="amber99sb-ildn",
                water_model="tip3p",
                ligand_resname="LIG",
                complex_gro=out / "complex.gro",
                complex_pdb=out / "complex.pdb",
                ligand_gro=out / "LIG_GMX.gro",
                ligand_itp=out / "LIG_GMX.itp",
                ligand_posre=out / "posre_LIG.itp",
                topol_top=out / "topol.top",
            )
            backend = GromacsBackend("native", ("gmx",), "gmx")
            simulator = BiophysicalMDSimulator(SimulationConfig(run_id="index_fail"), out, backend=backend)
            simulator._gmx = MagicMock(return_value=(0, "", ""))
            simulator._mdrun = MagicMock(return_value=None)

            with patch("md.simulation.write_index_from_gro", side_effect=ValueError("ligand missing")):
                with self.assertRaises(SimulationError):
                    simulator.run_pipeline(topology)

            manifest = json.loads((out / "md_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "failed")
            self.assertEqual(manifest["failed_stage"], "index_generation")
            self.assertIn("ligand missing", manifest["error"])


class TestProtocolBundleCompleteness(unittest.TestCase):
    def test_export_bundle_contains_replicate_mdps_and_loops(self):
        from md.engine import generate_standalone_scripts
        from md.simulation import SimulationConfig
        from md.ui import _build_mdp_bundle

        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            cfg = SimulationConfig(n_replicates=2)
            _build_mdp_bundle(out, cfg)
            self.assertTrue((out / "nvt_rep2.mdp").is_file())
            self.assertIn(
                "nsteps                  = 50000",
                (out / "nvt_rep2.mdp").read_text(encoding="utf-8"),
            )
            scripts = generate_standalone_scripts(out, n_replicates=2)
            self.assertIn("nvt_rep$i.mdp", scripts["bash"].read_text(encoding="utf-8"))
            self.assertIn("rep%%i_prod_md.tpr", scripts["batch"].read_text(encoding="utf-8"))

    def test_export_bundle_contains_ions_and_real_index_generator(self):
        from md.engine import generate_standalone_scripts
        from md.simulation import SimulationConfig
        from md.ui import _build_mdp_bundle

        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            cfg = SimulationConfig(water_model="tip4pew")
            _build_mdp_bundle(out, cfg)
            for name in ("ions.mdp", "em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"):
                self.assertTrue((out / name).is_file(), name)

            scripts = generate_standalone_scripts(out, water_model="tip4pew")
            bash = scripts["bash"].read_text(encoding="utf-8")
            self.assertIn("-cs tip4p.gro", bash)
            self.assertIn("make_index.py solv_ions.gro", bash)
            self.assertIn("trjconv", bash)


class TestAuditRound2Regressions(unittest.TestCase):
    def test_interaction_keys_and_chain_aware_mutation_engagement(self):
        from docking.interactions import compare_interactions

        healthy = {
            "polar_contact_candidates": [
                {"residue": "GLN102", "res_num": "102", "chain": "A",
                 "rec_atom": "OE1", "lig_atom": "OAB", "distance_angstrom": 2.9}
            ],
            "interacting_residues": [
                {"res_name": "GLN", "res_num": "102", "chain": "A", "min_distance": 2.9}
            ],
            "polar_contact_count": 1,
        }
        mutant = {"polar_contact_candidates": [], "interacting_residues": [], "polar_contact_count": 0}
        result = compare_interactions(healthy, mutant, ["A:102"])
        self.assertEqual(result["healthy_polar_contact_count"], 1)
        self.assertEqual(result["lost_polar_contact_residues"], ["A:GLN102"])
        self.assertTrue(result["mutation_engaged_in_contacts"])

    def test_pdbqt_aromatic_carbon_is_mapped_to_carbon(self):
        from docking.interactions import _parse_atoms_from_block

        line = "ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00    +0.000 A "
        atoms = _parse_atoms_from_block(line, is_ligand=True)
        self.assertEqual(atoms[0]["element"], "C")

    def test_timeseries_summary_has_sem_and_block_sem(self):
        from md.analysis import summarize_timeseries

        points = [(float(i), float(i)) for i in range(1, 21)]
        summary = summarize_timeseries(points, burn_in_fraction=0.0, n_blocks=4)
        self.assertEqual(int(summary["n"]), 20)
        self.assertAlmostEqual(summary["mean"], 10.5, places=3)
        self.assertGreater(summary["sd"], 0.0)
        self.assertGreater(summary["sem"], 0.0)
        self.assertGreaterEqual(summary["block_sem"], 0.0)

    def test_rcsb_quality_metadata_parser(self):
        from unittest.mock import MagicMock

        from docking.rcsb import fetch_entry_quality

        response = MagicMock()
        response.status_code = 200
        response.json.return_value = {
            "struct": {"title": "test"},
            "rcsb_entry_info": {"resolution_combined": [3.1], "experimental_method": "X-ray"},
            "refine": [{"ls_R_factor_R_free": 0.249, "ls_R_factor_R_work": 0.207}],
        }
        with patch("docking.rcsb.requests.get", return_value=response):
            quality = fetch_entry_quality("2JIT")
        self.assertEqual(quality["resolution_angstrom"], 3.1)
        self.assertEqual(quality["r_free"], 0.249)
        self.assertFalse(quality["meets_resolution_target"])
        self.assertTrue(quality["meets_rfree_target"])

    def test_drug_properties_include_pains_flag(self):
        from docking.interactions import calculate_drug_properties

        props = calculate_drug_properties("CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertIn("pains_alert", props)
        self.assertIn("pains_name", props)
        self.assertIn("aggregation_risk", props)

    def test_insertion_codes_do_not_collide_in_sasa_or_sequence_maps(self):
        from analysis.structure import calculate_sasa_map, get_protein_sequence, process_protein_structure

        lines = []
        serial = 1
        for resnum, icode in ((100, " "), (100, "A")):
            for name, x in (("N", 0.0), ("CA", 1.45), ("C", 2.0), ("O", 3.0)):
                lines.append(
                    f"ATOM  {serial:5d} {name:>4s} ALA A{resnum:4d}{icode:1s}   "
                    f"{x:8.3f}{float(serial):8.3f}{0.0:8.3f}  1.00 20.00           "
                    f"{name[0]:>2s}"
                )
                serial += 1
        lines.append("END")
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = Path(tmpdir) / "insertion.pdb"
            pdb_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            structure = process_protein_structure(pdb_path.read_text(encoding="utf-8"), "insertion")
        sequence = get_protein_sequence(structure, "A")
        self.assertEqual([residue["res_id"] for residue in sequence], ["100", "100A"])
        self.assertEqual(set(calculate_sasa_map(structure, "A").keys()), {"100", "100A"})

    def test_solute_diameter_is_derived_from_coordinates(self):
        from md.simulation import solute_bounding_diameter_nm

        gro = (
            "complex\n"
            "    2\n"
            + f"{1:5d}{'ALA':<5s}{'CA':>5s}{1:5d}{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}\n"
            + f"{1:5d}{'ALA':<5s}{'CA':>5s}{2:5d}{2.5:8.3f}{1.0:8.3f}{0.5:8.3f}\n"
            + "   5.0 5.0 5.0\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "complex.gro"
            path.write_text(gro, encoding="utf-8")
            diameter = solute_bounding_diameter_nm(path)
        self.assertAlmostEqual(diameter, 2.5, places=2)

    def test_replicate_nvt_uses_configured_duration_not_misordered_args(self):
        from unittest.mock import MagicMock

        from md.engine import GromacsBackend
        from md.simulation import BiophysicalMDSimulator, SimulationConfig
        from md.topology import TopologyResult

        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            for name in ("complex.gro", "complex.pdb", "LIG_GMX.gro", "LIG_GMX.itp", "posre_LIG.itp",
                         "npt.mdp", "md.mdp"):
                (out / name).write_text("placeholder\n", encoding="utf-8")
            topology = TopologyResult(
                method="test", force_field="amber99sb-ildn", water_model="tip3p",
                ligand_resname="LIG", complex_gro=out / "complex.gro",
                complex_pdb=out / "complex.pdb", ligand_gro=out / "LIG_GMX.gro",
                ligand_itp=out / "LIG_GMX.itp", ligand_posre=out / "posre_LIG.itp",
                topol_top=out / "topol.top",
            )
            cfg = SimulationConfig(nvt_duration_ps=0.005, npt_duration_ps=0.005,
                                   prod_duration_ns=0.0005, n_replicates=2)
            simulator = BiophysicalMDSimulator(
                cfg, out, backend=GromacsBackend("native", ("gmx",), "gmx")
            )
            simulator._gmx = MagicMock(return_value=(0, "", ""))
            simulator._mdrun = MagicMock(return_value=None)
            simulator._replicate_summary = MagicMock(return_value={"replicate": 2})
            simulator._run_extra_replicate(
                2, topology, {"npt.mdp": out / "npt.mdp", "md.mdp": out / "md.mdp"}
            )
            nvt_text = (out / "rep2_nvt.mdp").read_text(encoding="utf-8")
        self.assertIn("nsteps                  = 2", nvt_text)
        self.assertIn("gen_seed                = 2042", nvt_text)

    def test_mmpbsa_parser_reads_modern_delta_output(self):
        from md.analysis import load_mmpbsa_results

        sample = (
            "gmx_MMPBSA Version=v1.7.0\n"
            "Calculations performed using 16 complex frames\n"
            "igb = 5, intdiel = 2.0\n"
            "Complex:\n"
            "VDWAALS -2035.53\n"
            "EEL -16750.43\n"
            "Delta (Complex - Receptor - Ligand):\n"
            "ΔVDWAALS -58.89\n"
            "ΔEEL -31.13\n"
            "ΔEGB 40.99\n"
            "ΔESURF -8.22\n"
            "ΔTOTAL -57.25\n"
            "ΔG binding =  -48.0704 +/-  2.1972\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "FINAL_RESULTS_MMPBSA.dat").write_text(sample, encoding="utf-8")
            result = load_mmpbsa_results(Path(tmpdir))
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result.delta_g_bind, -48.07, places=2)
        self.assertAlmostEqual(result.delta_vdw, -58.89, places=2)
        self.assertAlmostEqual(result.delta_elec, -31.13, places=2)
        self.assertEqual(result.frames_evaluated, 16)
        self.assertEqual(result.igb, 5)
        self.assertEqual(result.solvation_model, "MM-GBSA (OBC2)")
        self.assertTrue(result.internal_dielectric_specified)
        self.assertFalse(result.entropy_omitted)

    def test_hetatm_mse_is_converted_not_dropped(self):
        from md.topology import clean_receptor_for_md

        pdb = (
            "HETATM    1  N   MSE A   1       0.000   0.000   0.000  1.00 20.00           N\n"
            "HETATM    2  SE  MSE A   1       1.000   0.000   0.000  1.00 20.00          SE\n"
            "ATOM      3  CA  ALA A   2       2.000   0.000   0.000  1.00 20.00           C\n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            src = Path(tmpdir) / "in.pdb"
            src.write_text(pdb, encoding="utf-8")
            output, _ = clean_receptor_for_md(src, Path(tmpdir) / "out.pdb")
            text = output.read_text(encoding="utf-8")
        self.assertIn("MET", text)
        self.assertIn(" SD ", text)
        self.assertNotIn("MSE", text)

    def test_contact_map_is_protein_only_by_default(self):
        from analysis.structure import calculate_all_distances, process_protein_structure

        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 20.00           C\n"
            "ATOM      4  O   ALA A   1       3.000   1.400   0.000  1.00 20.00           O\n"
            "ATOM      5  N   ALA A  10       0.000  10.000   0.000  1.00 20.00           N\n"
            "ATOM      6  CA  ALA A  10       1.450  10.000   0.000  1.00 20.00           C\n"
            "ATOM      7  C   ALA A  10       2.000  11.400   0.000  1.00 20.00           C\n"
            "ATOM      8  O   ALA A  10       3.000  11.400   0.000  1.00 20.00           O\n"
            "HETATM    9  FE  HEM A 101       1.000   0.500   0.500  1.00 20.00          FE\n"
            "END\n"
        )
        structure = process_protein_structure(pdb, "contacts_test")
        protein_only = calculate_all_distances(structure, "A", radius=5.0)
        included = calculate_all_distances(structure, "A", radius=5.0, include_hetero=True)
        self.assertEqual(protein_only[0]["min_dist"], "-")
        self.assertNotEqual(included[0]["min_dist"], "-")

    def test_aggregator_structural_alerts_are_transparent(self):
        from docking.interactions import calculate_drug_properties

        quinone = calculate_drug_properties("O=C1C=CC(=O)C=C1")
        self.assertTrue(quinone["aggregation_risk"])
        self.assertIn("quinone", quinone["aggregation_alerts"])
        self.assertEqual(quinone["aggregation_logp_threshold"], 3.5)
        ethanol = calculate_drug_properties("CCO")
        self.assertFalse(ethanol["aggregation_risk"])
        self.assertEqual(ethanol["aggregation_alerts"], [])

    def test_low_rscc_residues_are_listed_for_the_ui(self):
        from unittest.mock import MagicMock, patch

        from docking.rcsb import fetch_entry_quality

        xml = b"""<?xml version="1.0"?>
        <wwPDB-validation-information>
          <ModelledSubgroup chain="A" resnum="1" icode=" " resname="ALA" rscc="0.95"/>
          <ModelledSubgroup chain="A" resnum="2" icode=" " resname="GLY" rscc="0.70"/>
          <ModelledSubgroup chain="A" resnum="3" icode=" " resname="MET" rscc="0.85"/>
        </wwPDB-validation-information>"""
        response = MagicMock(status_code=200)
        response.json.return_value = {
            "rcsb_entry_info": {"resolution_combined": [1.8], "experimental_method": "X-RAY"},
            "refine": [{"ls_R_factor_R_free": 0.2, "ls_R_factor_R_work": 0.18}],
            "struct": {"title": "test"},
        }
        with patch("docking.rcsb.requests.get", return_value=response), \
             patch("docking.rcsb.fetch_validation_xml", return_value=xml):
            quality = fetch_entry_quality("3HTB", include_residue_rscc=True)
        self.assertEqual(quality["residue_rscc_summary"]["n_residues"], 3)
        self.assertEqual(quality["residue_rscc_low"][0]["resnum"], "2")
        self.assertEqual(quality["residue_rscc_low"][0]["rscc"], 0.7)

    def test_residue_rscc_summary_is_conservative(self):
        from docking.rcsb import parse_residue_rscc, summarize_residue_rscc

        xml = b"""<?xml version="1.0"?>
        <wwPDB-validation-information>
          <ModelledSubgroup chain="A" resnum="1" icode=" " resname="ALA" rscc="0.95"/>
          <ModelledSubgroup chain="A" resnum="2" icode=" " resname="GLY" rscc="0.75"/>
          <ModelledSubgroup chain="A" resnum="2" icode=" " resname="GLY" rscc="0.70"/>
          <ModelledSubgroup chain="A" resnum="3" icode=" " resname="MET" rscc="0.85"/>
          <ModelledSubgroup chain="A" resnum="101" icode=" " resname="JZ4" rscc="0.95"/>
        </wwPDB-validation-information>"""
        residues = parse_residue_rscc(xml)
        summary = summarize_residue_rscc(residues)
        self.assertEqual(summary["n_residues"], 3)
        self.assertEqual(summary["min_rscc"], 0.70)
        self.assertEqual(summary["median_rscc"], 0.85)
        self.assertEqual(summary["n_below_0_8"], 1)

    def test_ligand_rscc_parser_uses_validation_report(self):
        from docking.rcsb import parse_ligand_rscc

        xml = b"""<?xml version="1.0"?>
        <wwPDB-validation-information>
          <ModelledSubgroup chain="A" resnum="167" icode=" " resname="JZ4"
                            rscc="0.953" avgoccu="1.000"/>
          <ModelledSubgroup chain="B" resnum="167" icode=" " resname="JZ4"
                            rscc="0.812" avgoccu="0.900"/>
        </wwPDB-validation-information>"""
        parsed = parse_ligand_rscc(xml, "jz4")
        self.assertEqual(parsed["ligand_code"], "JZ4")
        self.assertEqual(parsed["ligand_rscc"], 0.812)
        self.assertEqual(len(parsed["entries"]), 2)

    def test_tip4p_ew_water_topology_is_validated(self):
        from md.engine import GromacsBackend
        from md.topology import inspect_water_topology_text, validate_water_model_topology

        itp = (
            "[ moleculetype ]\nSOL 2\n[ atoms ]\n"
            "1 OW_tip4pew 1 SOL OW 1 0 16.0\n"
            "2 HW_tip4pew 1 SOL HW1 1 0.52422 1.008\n"
            "3 HW_tip4pew 1 SOL HW2 1 0.52422 1.008\n"
            "4 MW 1 SOL MW 1 -1.04844 0.0\n"
        )
        gro = "tip4p\n 864\n" + "".join(
            f"{1:5d}{'SOL':<5s}{name:<5s}{i:5d}{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}\n"
            for i, name in enumerate(["OW", "HW1", "HW2", "MW"], start=1)
        ) + "1.8 1.8 1.8\n"
        self.assertTrue(inspect_water_topology_text(itp, gro, "tip4pew")["validated"])
        self.assertFalse(
            inspect_water_topology_text(itp.replace("4 MW", ";4 MW"), gro, "tip4pew")["validated"]
        )

        backend = GromacsBackend("wsl", ("wsl", "-e", "/usr/bin/gmx"), "/usr/bin/gmx", data_prefix="/usr")
        with patch("md.topology.run_backend_shell", return_value=(0, itp + "@@@WATER@@@" + gro, "")):
            result = validate_water_model_topology(backend, "amber99sb-ildn", "tip4pew")
        self.assertTrue(result["validated"])
        self.assertEqual(result["moleculetype"], "SOL")

    def test_crystal_water_retention_builds_complete_water(self):
        import math

        from md.indexing import parse_gro_atoms
        from md.topology import (
            append_crystal_waters_to_gro,
            append_crystal_waters_to_topology,
            extract_crystal_waters,
        )

        pdb = (
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
            "HETATM    2  O   HOH A 501       0.500   0.500   0.500  1.00 20.00           O\n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            src = base / "water.pdb"
            src.write_text(pdb, encoding="utf-8")
            waters = extract_crystal_waters(src)
            self.assertEqual(len(waters), 1)
            water = waters[0]
            self.assertEqual(water["resname"], "SOL")
            for key in ("ow", "hw1", "hw2"):
                self.assertIn(key, water)
            for h in ("hw1", "hw2"):
                distance = math.dist(water["ow"], water[h])
                self.assertAlmostEqual(distance, 0.09572, places=4)
            angle_vec_1 = [water["hw1"][i] - water["ow"][i] for i in range(3)]
            angle_vec_2 = [water["hw2"][i] - water["ow"][i] for i in range(3)]
            dot = sum(a * b for a, b in zip(angle_vec_1, angle_vec_2))
            norm = math.dist(water["hw1"], water["ow"]) * math.dist(water["hw2"], water["ow"])
            self.assertAlmostEqual(math.degrees(math.acos(dot / norm)), 104.52, places=1)

            top = base / "topol.top"
            top.write_text("[ molecules ]\nProtein 1\nLIG 1\n", encoding="utf-8")
            append_crystal_waters_to_topology(top, len(waters))
            self.assertIn("SOL", top.read_text(encoding="utf-8"))

            gro = base / "complex.gro"
            gro.write_text(
                "complex\n    1\n"
                f"{1:5d}{'ALA':<5s}{'CA':<5s}{1:5d}{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}\n"
                "   1.0 1.0 1.0\n",
                encoding="utf-8",
            )
            append_crystal_waters_to_gro(gro, waters)
            atoms, _ = parse_gro_atoms(gro)
        self.assertEqual(len(atoms), 4)
        self.assertEqual([atom[2] for atom in atoms[-3:]], ["OW", "HW1", "HW2"])

    def test_monatomic_ion_retention_helpers(self):
        from md.indexing import parse_gro_atoms
        from md.topology import (
            append_ions_to_gro,
            append_ions_to_topology,
            extract_monatomic_ions,
        )

        pdb = (
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
            "HETATM    2 ZN    ZN A 201       1.000   2.000   3.000  1.00 20.00          ZN\n"
            "HETATM    3  FE  HEM A 301       4.000   5.000   6.000  1.00 20.00          FE\n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            src = base / "receptor.pdb"
            src.write_text(pdb, encoding="utf-8")
            ions, unsupported = extract_monatomic_ions(src)
            self.assertEqual(len(ions), 1)
            self.assertEqual(ions[0]["resname"], "ZN")
            self.assertEqual(unsupported, {"HEM"})

            top = base / "topol.top"
            top.write_text(
                "[ molecules ]\nProtein_chain_A 1\nLIG 1\n", encoding="utf-8"
            )
            append_ions_to_topology(top, ions)
            self.assertIn("ZN                   1", top.read_text(encoding="utf-8"))

            gro = base / "complex.gro"
            gro.write_text(
                "complex\n    1\n"
                f"{1:5d}{'ALA':<5s}{'CA':<5s}{1:5d}{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}\n"
                "   1.0 1.0 1.0\n",
                encoding="utf-8",
            )
            append_ions_to_gro(gro, ions)
            atoms, _ = parse_gro_atoms(gro)
        self.assertEqual(len(atoms), 2)
        self.assertEqual(atoms[-1][1], "ZN")

    def test_sasa_is_protein_only_but_hetero_is_retained_for_display(self):
        from analysis.structure import calculate_sasa_map, process_protein_structure

        pdb = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 20.00           C\n"
            "ATOM      4  O   ALA A   1       3.000   1.400   0.000  1.00 20.00           O\n"
            "HETATM    5  FE  HEM A 101       1.500   0.500   0.500  1.00 20.00          FE\n"
            "END\n"
        )
        structure = process_protein_structure(pdb, "sasa_test")
        chain = structure[0]["A"]
        self.assertIn("HEM", [residue.get_resname() for residue in chain])
        self.assertEqual(set(calculate_sasa_map(structure, "A").keys()), {"1"})
        self.assertTrue(hasattr(next(iter(chain)), "sasa"))

    def test_critical_cofactor_removal_is_fail_closed(self):
        from md.topology import TopologyBuildError, clean_receptor_for_md

        pdb = (
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
            "HETATM    2  FE  HEM A 101       1.000   0.000   0.000  1.00 20.00          FE\n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            src = Path(tmpdir) / "heme.pdb"
            src.write_text(pdb, encoding="utf-8")
            with self.assertRaises(TopologyBuildError):
                clean_receptor_for_md(src, Path(tmpdir) / "out.pdb")
            _, warnings = clean_receptor_for_md(
                src, Path(tmpdir) / "out.pdb", allow_cofactor_removal=True
            )
        self.assertTrue(any("HEM" in warning for warning in warnings))

    def test_receptor_sanitizer_selects_first_model_and_best_altloc(self):
        from docking.receptor_prep import prepare_receptor

        pdb = (
            "MODEL        1\n"
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  0.50 20.00           N\n"
            "ATOM      2  CA AALA A   1       1.450   0.000   0.000  0.50 20.00           C\n"
            "ATOM      3  CA BALA A   1       9.000   0.000   0.000  0.80 20.00           C\n"
            "ENDMDL\n"
            "MODEL        2\n"
            "ATOM      1  N   ALA A   1       9.000   0.000   0.000  1.00 20.00           N\n"
            "ENDMDL\nEND\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            clean, _, report = prepare_receptor(pdb, Path(tmpdir))
            text = clean.read_text(encoding="utf-8")
        self.assertEqual(report["models_dropped"], 1)
        self.assertEqual(report["altloc_records_resolved"], 2)
        self.assertNotIn("MODEL", text)
        ca_line = next(line for line in text.splitlines() if line[12:16].strip() == "CA")
        self.assertAlmostEqual(float(ca_line[30:38]), 9.0, places=2)

    def test_native_windows_force_field_discovery(self):
        from md.engine import GromacsBackend, list_force_fields, resolve_force_field

        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "share" / "gromacs" / "top" / "amber99sb-ildn.ff").mkdir(parents=True)
            backend = GromacsBackend("native", ("gmx",), "gmx", data_prefix=tmpdir)
            self.assertIn("amber99sb-ildn", list_force_fields(backend))
            self.assertEqual(resolve_force_field(backend, "amber99sb-ildn"), "amber99sb-ildn")

    def test_matched_comparison_requires_pair_validation(self):
        from docking.comparison import run_matched_docking_comparison

        with self.assertRaisesRegex(ValueError, "validated homologous"):
            run_matched_docking_comparison(
                healthy_pdb_text="ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n",
                mutant_pdb_text="ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n",
                healthy_target_meta={"pdb_id": "WT"},
                mutant_target_meta={"pdb_id": "MUT"},
                ligand_smiles="CC",
                ligand_name="Ethanol",
                healthy_pocket_residues=["A:1"],
                mutant_pocket_residues=["A:1"],
            )

    def test_empty_alignment_returns_invalid_metrics(self):
        metrics = calculate_sequence_identity("", "")
        self.assertFalse(metrics["valid"])
        self.assertEqual(metrics["identity_pct"], 0.0)


if __name__ == "__main__":
    unittest.main()
