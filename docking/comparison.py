# Matched Healthy-versus-Mutant Docking Comparison Module
import datetime
import hashlib
from io import StringIO
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D

from docking.models import GridBox, PoseCluster
from docking.engine import get_vina_version, run_vina_multi_seeds
from docking.binding_site import define_grid_from_residues
from docking.receptor_prep import prepare_receptor
from docking.ligand_prep import prepare_ligand_from_smiles
from docking.validation import (
    cluster_poses_across_seeds,
    validate_redocking_top_n,
    verify_reference_ligand_identity,
    pose_block_to_mol,
)
from docking.rcsb import fetch_chem_comp_smiles
from docking.storage import (
    create_run_directory,
    write_initial_manifest,
    build_and_save_manifest,
    compute_file_sha256,
    DISCLAIMER_TEXT,
)


COMPARISON_LABEL = (
    "Protocol-specific docking-score ranking difference. "
    "It is not a measured affinity difference, drug-response prediction, or clinical conclusion."
)

PROTONATION_ASSUMPTION_NOTE = (
    "The submitted SMILES is protonated at pH 7.4 with OpenBabel's pH-aware model "
    "before docking; tautomers are not enumerated. The manifest records the exact "
    "prepared protomer, formal charge and method."
)

MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT: int = 2


def validate_pocket_chain(residue_specs: list[str], selected_chain: str, label: str) -> None:
    for spec in residue_specs:
        if ":" not in spec:
            raise ValueError(
                f"{label} pocket residue {spec!r} must use CHAIN:RESIDUE notation."
            )
        chain, _ = spec.split(":", 1)
        if chain.strip() != selected_chain:
            raise ValueError(
                f"{label} pocket residue {spec!r} does not belong to selected chain {selected_chain!r}."
            )

def _parse_pocket_spec(spec: str) -> tuple[str | None, int, str]:
    """Parse ``CHAIN:RESIDUE[ICODE]`` into (chain, resseq, icode)."""
    raw = str(spec).strip()
    if ":" in raw:
        chain, residue = (part.strip() for part in raw.split(":", 1))
    else:
        chain, residue = None, raw
    digits = "".join(ch for ch in residue if ch.isdigit())
    icode = "".join(ch for ch in residue if ch.isalpha())
    if not digits:
        raise ValueError(f"Invalid pocket residue specification: {spec!r}")
    return chain, int(digits), icode.upper()


def _protein_ca_map(pdb_path: Path) -> dict[tuple[str, int, str], np.ndarray]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_path))
    model = structure[0]
    coords: dict[tuple[str, int, str], np.ndarray] = {}
    for chain in model:
        chain_id = str(chain.id)
        for residue in chain:
            if residue.id[0] != " " or "CA" not in residue:
                continue
            key = (chain_id, int(residue.id[1]), str(residue.id[2]).strip().upper())
            coords[key] = np.asarray(residue["CA"].coord, dtype=float)
    return coords


def superimpose_receptors(
    healthy_pdb: Path,
    mutant_pdb: Path,
    min_common_ca: int = 30,
) -> dict | None:
    """Superimpose two receptor structures on shared C-alpha atoms.

    Matched docking compares a WT and a variant of the same protein, so a rigid
    C-alpha superposition is valid.  Returns ``None`` when fewer than
    ``min_common_ca`` equivalent residues exist (e.g. different chain naming).
    """
    try:
        healthy_ca = _protein_ca_map(healthy_pdb)
        mutant_ca = _protein_ca_map(mutant_pdb)
    except Exception:
        return None
    common = sorted(set(healthy_ca) & set(mutant_ca))
    if len(common) < min_common_ca:
        return None

    fixed = np.asarray([healthy_ca[key] for key in common], dtype=float)
    moving = np.asarray([mutant_ca[key] for key in common], dtype=float)
    try:
        superimposer = SVDSuperimposer()
        superimposer.set(fixed, moving)
        superimposer.run()
        rotation, translation = superimposer.get_rotran()
        rmsd = float(superimposer.get_rms())
    except Exception:
        return None
    return {
        "rotation": np.asarray(rotation, dtype=float),
        "translation": np.asarray(translation, dtype=float),
        "receptor_rmsd_angstrom": round(rmsd, 3),
        "common_ca": len(common),
    }


def _pocket_centroid(pdb_path: Path, residue_specs: list[str] | None) -> np.ndarray | None:
    if not residue_specs:
        return None
    try:
        parser = PDBParser(QUIET=True)
        model = parser.get_structure("receptor", str(pdb_path))[0]
        wanted = {_parse_pocket_spec(spec) for spec in residue_specs}
        points = []
        for chain in model:
            for residue in chain:
                key = (str(chain.id), int(residue.id[1]), str(residue.id[2]).strip().upper())
                if key in wanted:
                    points.extend(atom.coord for atom in residue.get_atoms())
        if not points:
            return None
        return np.asarray(points, dtype=float).mean(axis=0)
    except Exception:
        return None


def _pose_centroid(pose_block: str) -> np.ndarray | None:
    mol = pose_block_to_mol(pose_block)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    conformer = mol.GetConformer()
    points = [
        [conformer.GetAtomPosition(i).x, conformer.GetAtomPosition(i).y, conformer.GetAtomPosition(i).z]
        for i in range(mol.GetNumAtoms())
    ]
    return np.asarray(points, dtype=float).mean(axis=0)


def evaluate_binding_mode_compatibility(
    healthy_pose_block: str,
    mutant_pose_block: str,
    healthy_receptor_pdb: Path | None = None,
    mutant_receptor_pdb: Path | None = None,
    healthy_pocket_residues: list[str] | None = None,
    mutant_pocket_residues: list[str] | None = None,
    rmsd_cutoff: float = 2.0,
) -> tuple[bool | None, dict]:
    """Check that both consensus poses describe the same binding mode.

    Returns ``(compatible, details)`` where ``compatible`` is ``True``,
    ``False`` (a real incompatibility), or ``None`` when the check could not be
    evaluated (unparseable pose).  A ``None`` result is reported as a warning
    rather than treated as evidence of incompatibility.

    Primary test: superpose the mutant receptor onto the healthy receptor and
    compare the mutant ligand pose with the healthy ligand pose using RDKit's
    symmetry-aware, no-fit ``CalcRMS``.  When a receptor superposition is not
    possible, fall back to a weaker geometric check (ligand centroid must sit at
    the mapped pocket centroid) and record that the stronger test was skipped.
    """
    details: dict = {
        "rmsd_cutoff_angstrom": rmsd_cutoff,
        "method": "receptor_CA_superposition + symmetry-aware ligand RMSD",
    }

    healthy_pose_mol = pose_block_to_mol(healthy_pose_block)
    mutant_pose_mol = pose_block_to_mol(mutant_pose_block)
    if healthy_pose_mol is None or mutant_pose_mol is None:
        details["reason"] = "Could not parse one or both ligand poses."
        details["evaluated"] = False
        return None, details
    if healthy_pose_mol.GetNumAtoms() != mutant_pose_mol.GetNumAtoms():
        details["reason"] = (
            f"Ligand heavy-atom counts differ ({healthy_pose_mol.GetNumAtoms()} vs "
            f"{mutant_pose_mol.GetNumAtoms()})."
        )
        details["evaluated"] = True
        return False, details

    superposition = None
    if healthy_receptor_pdb and mutant_receptor_pdb:
        superposition = superimpose_receptors(Path(healthy_receptor_pdb), Path(mutant_receptor_pdb))

    if superposition is not None:
        rotation = superposition["rotation"]
        translation = superposition["translation"]
        transformed = Chem.Mol(mutant_pose_mol)
        conformer = transformed.GetConformer()
        for index in range(transformed.GetNumAtoms()):
            point = conformer.GetAtomPosition(index)
            coord = np.dot(np.asarray([point.x, point.y, point.z]), rotation) + translation
            conformer.SetAtomPosition(index, Point3D(float(coord[0]), float(coord[1]), float(coord[2])))
        try:
            pose_rmsd = float(rdMolAlign.CalcRMS(transformed, healthy_pose_mol))
        except Exception as exc:
            details["reason"] = f"Symmetry-aware ligand RMSD failed: {exc}"
            return False, details
        details.update({
            "receptor_superposition_rmsd_angstrom": superposition["receptor_rmsd_angstrom"],
            "common_ca": superposition["common_ca"],
            "ligand_pose_rmsd_angstrom": round(pose_rmsd, 3),
            "evaluated": True,
            "method": "receptor_CA_superposition + symmetry-aware ligand RMSD",
        })
        return pose_rmsd <= rmsd_cutoff, details

    # Fallback: the two poses must at least occupy the homologous pocket.
    details["method"] = "pocket_centroid_geometric_fallback (no receptor superposition)"
    healthy_centroid = _pocket_centroid(Path(healthy_receptor_pdb), healthy_pocket_residues) if healthy_receptor_pdb else None
    mutant_centroid = _pocket_centroid(Path(mutant_receptor_pdb), mutant_pocket_residues) if mutant_receptor_pdb else None
    healthy_ligand_centroid = _pose_centroid(healthy_pose_block)
    mutant_ligand_centroid = _pose_centroid(mutant_pose_block)
    if healthy_centroid is None or mutant_centroid is None or healthy_ligand_centroid is None or mutant_ligand_centroid is None:
        details["reason"] = "Neither receptor superposition nor pocket centroids were available."
        details["evaluated"] = False
        return None, details
    healthy_offset = float(np.linalg.norm(healthy_ligand_centroid - healthy_centroid))
    mutant_offset = float(np.linalg.norm(mutant_ligand_centroid - mutant_centroid))
    details.update({
        "healthy_ligand_to_pocket_centroid_angstrom": round(healthy_offset, 3),
        "mutant_ligand_to_pocket_centroid_angstrom": round(mutant_offset, 3),
        "pocket_centroid_offset_limit_angstrom": 10.0,
        "evaluated": True,
    })
    compatible = healthy_offset <= 10.0 and mutant_offset <= 10.0
    if not compatible:
        details["reason"] = "At least one ligand pose lies outside the mapped binding pocket."
    return compatible, details


def validate_matched_conditions(
    healthy_domain: str,
    mutant_domain: str,
    healthy_ligand_hash: str,
    mutant_ligand_hash: str,
    healthy_waters_policy: str,
    mutant_waters_policy: str,
    healthy_seeds: list,
    mutant_seeds: list,
    healthy_exhaustiveness: int,
    mutant_exhaustiveness: int
) -> tuple[bool, str]:
    """Strictly verify matched comparison requirements before allowing delta calculations."""
    if healthy_domain != mutant_domain:
        return False, f"Domain mismatch: healthy domain '{healthy_domain}' != mutant domain '{mutant_domain}'."
    if healthy_ligand_hash != mutant_ligand_hash:
        return False, "Ligand mismatch: prepared ligand SHA-256 hashes must be identical."
    if healthy_waters_policy != mutant_waters_policy:
        return False, "Receptor preparation mismatch: waters policy differs."
    if healthy_seeds != mutant_seeds or len(healthy_seeds) < 3:
        return False, "Seed protocol mismatch: must run identical 3 independent seeds."
    if healthy_exhaustiveness != mutant_exhaustiveness:
        return False, f"Search exhaustiveness mismatch: {healthy_exhaustiveness} != {mutant_exhaustiveness}."
    return True, "Matched comparison conditions satisfied."


def _select_consensus_cluster(clusters: list[PoseCluster]) -> PoseCluster:
    if not clusters:
        raise ValueError("Matched comparison rejected: no pose clusters were produced.")
    # Prefer reproducibility first, then a lower median Vina score.
    return sorted(
        clusters,
        key=lambda cluster: (-cluster.seed_count, cluster.median_score, cluster.top_score),
    )[0]


def verify_comparison_clusters(
    h_clusters: list[PoseCluster],
    m_clusters: list[PoseCluster],
    healthy_receptor_pdb: Path | None = None,
    mutant_receptor_pdb: Path | None = None,
    healthy_pocket_residues: list[str] | None = None,
    mutant_pocket_residues: list[str] | None = None,
    binding_mode_cutoff: float = 2.0,
    details: dict | None = None,
) -> tuple[PoseCluster, PoseCluster, float, list[str]]:
    """Select the consensus cluster of each arm and require a shared binding mode.

    A ΔScore is only meaningful when both consensus poses occupy the same
    binding mode.  The gate compares the two representative poses after
    superposing the mutant receptor onto the healthy receptor (symmetry-aware,
    no-fit ligand RMSD) and refuses the comparison otherwise.
    """
    healthy_cluster = _select_consensus_cluster(h_clusters)
    mutant_cluster = _select_consensus_cluster(m_clusters)

    if healthy_cluster.seed_count < 2 or mutant_cluster.seed_count < 2:
        raise ValueError(
            "Matched comparison rejected: the selected consensus pose must be supported "
            "by at least two independent seeds for both receptors. "
            f"healthy={healthy_cluster.seed_count}, mutant={mutant_cluster.seed_count}."
        )

    warnings: list[str] = []
    if details is not None:
        details.clear()
    if healthy_receptor_pdb and mutant_receptor_pdb:
        compatible, mode_details = evaluate_binding_mode_compatibility(
            healthy_cluster.representative_pose.pdbqt_block,
            mutant_cluster.representative_pose.pdbqt_block,
            healthy_receptor_pdb=healthy_receptor_pdb,
            mutant_receptor_pdb=mutant_receptor_pdb,
            healthy_pocket_residues=healthy_pocket_residues,
            mutant_pocket_residues=mutant_pocket_residues,
            rmsd_cutoff=binding_mode_cutoff,
        )
        if details is not None:
            details.update(mode_details)
        if compatible is False:
            raise ValueError(
                "Matched comparison rejected: the healthy and mutant consensus poses do not "
                "represent the same binding mode, so a ΔScore would be uninterpretable. "
                f"Details: {mode_details.get('reason', mode_details)}"
            )
        if compatible is None:
            message = (
                "Binding-mode compatibility could not be evaluated for this pair "
                f"({mode_details.get('reason', 'unknown reason')}); ΔScore must be treated as exploratory."
            )
            warnings.append(message)
            if details is not None:
                details["reason"] = message
    else:
        message = (
            "Binding-mode compatibility gate could not be evaluated because receptor "
            "coordinates were not supplied; ΔScore must be treated as exploratory."
        )
        warnings.append(message)
        if details is not None:
            details["reason"] = message
            details["method"] = "not_evaluated"

    delta_score = round(mutant_cluster.median_score - healthy_cluster.median_score, 3)
    return healthy_cluster, mutant_cluster, delta_score, warnings


def run_matched_docking_comparison(
    healthy_pdb_text: str,
    mutant_pdb_text: str,
    healthy_target_meta: dict,
    mutant_target_meta: dict,
    ligand_smiles: str,
    ligand_name: str,
    healthy_pocket_residues: list[str] | None = None,
    mutant_pocket_residues: list[str] | None = None,
    pocket_residues: list[str] | None = None,
    target_domain: str = "",
    healthy_domain: str = "",
    mutant_domain: str = "",
    padding: float = 8.0,
    seeds: list[int] | None = None,
    exhaustiveness: int = 8,
    input_fingerprint: str | None = None,
    pair_validation_meta: dict | None = None,
    wt_reference_ligand_meta: dict | None = None,
    pubchem_cid: int | None = None,
    keep_waters: bool = False,
) -> dict:
    """
    Execute controlled matched docking comparison between healthy and mutant structures.
    Refuses comparison unless all protocol settings match strictly and both dominant
    clusters have at least 2 independent supporting seeds.

    Accepts separate healthy and mutant pocket residue lists to accommodate chain ID
    differences or residue numbering shifts.
    """
    if seeds is None:
        seeds = [42, 101, 2024]
    if not healthy_pdb_text or not mutant_pdb_text:
        raise ValueError("Missing coordinates for healthy or mutant protein.")

    # Handle backward compatibility with single pocket_residues parameter
    if healthy_pocket_residues is None and pocket_residues is not None:
        healthy_pocket_residues = pocket_residues
    if mutant_pocket_residues is None and pocket_residues is not None:
        mutant_pocket_residues = pocket_residues

    if not healthy_pocket_residues or not mutant_pocket_residues:
        raise ValueError("Docking was not run: define validated healthy and mutant binding sites.")

    pair_meta = pair_validation_meta or {}
    if not pair_meta.get("is_valid"):
        raise ValueError(
            "Matched comparison rejected: a validated homologous WT/variant pair is required. "
            "Run validate_matched_pair() first and pass its metadata as pair_validation_meta."
        )

    if len(seeds) < 3:
        raise ValueError("Matched comparison requires at least 3 independent seeds.")

    vina_ver = get_vina_version()
    h_created_at = healthy_target_meta.get("created_at_utc") or datetime.datetime.now(datetime.timezone.utc).isoformat()
    m_created_at = mutant_target_meta.get("created_at_utc") or datetime.datetime.now(datetime.timezone.utc).isoformat()

    # 1. Setup independent run directories
    h_run_id, h_run_dir = create_run_directory("dock_matched_healthy")
    m_run_id, m_run_dir = create_run_directory("dock_matched_mutant")

    # Save immutable source PDB files
    with open(h_run_dir / "receptor_source.pdb", "w", encoding="utf-8") as f:
        f.write(healthy_pdb_text)
    with open(m_run_dir / "receptor_source.pdb", "w", encoding="utf-8") as f:
        f.write(mutant_pdb_text)

    # 2. Compute individual pocket grids using each receptor's own mapped residues
    h_grid = define_grid_from_residues(healthy_pdb_text, healthy_pocket_residues, padding=padding)
    m_grid = define_grid_from_residues(mutant_pdb_text, mutant_pocket_residues, padding=padding)

    grid_meta_healthy = {
        "center": [h_grid.center_x, h_grid.center_y, h_grid.center_z],
        "size": [h_grid.size_x, h_grid.size_y, h_grid.size_z],
    }
    grid_meta_mutant = {
        "center": [m_grid.center_x, m_grid.center_y, m_grid.center_z],
        "size": [m_grid.size_x, m_grid.size_y, m_grid.size_z],
    }

    engine_meta = {
        "name": "AutoDock Vina",
        "version": vina_ver,
        "seeds": seeds,
        "exhaustiveness": exhaustiveness,
        "num_modes": 9,
        "energy_range": 3.0,
    }
    binding_site_meta_healthy = {
        "method": "residues",
        "residues": healthy_pocket_residues,
        "padding": padding,
        "mapped_mutant_residues": mutant_pocket_residues,
        "pair_validation": pair_validation_meta or {},
    }
    binding_site_meta_mutant = {
        "method": "residues",
        "residues": mutant_pocket_residues,
        "padding": padding,
        "source_healthy_residues": healthy_pocket_residues,
        "pair_validation": pair_validation_meta or {},
    }

    ligand_meta_initial = {
        "name": ligand_name,
        "input_smiles": ligand_smiles,
        "pubchem_cid": pubchem_cid,
        "protonation_note": PROTONATION_ASSUMPTION_NOTE,
    }

    # 3. Create initial manifest for both runs BEFORE preparation or engine runs
    h_initial = {
        "schema_version": 1,
        "run_id": h_run_id,
        "status": "running",
        "created_at_utc": h_created_at,
        "target": dict(healthy_target_meta, pair_validation=pair_validation_meta or {}),
        "binding_site": binding_site_meta_healthy,
        "ligand": ligand_meta_initial,
        "grid": grid_meta_healthy,
        "engine": engine_meta,
    }
    write_initial_manifest(h_run_dir, h_initial)

    m_initial = {
        "schema_version": 1,
        "run_id": m_run_id,
        "status": "running",
        "created_at_utc": m_created_at,
        "target": dict(mutant_target_meta, pair_validation=pair_validation_meta or {}),
        "binding_site": binding_site_meta_mutant,
        "ligand": ligand_meta_initial,
        "grid": grid_meta_mutant,
        "engine": engine_meta,
    }
    write_initial_manifest(m_run_dir, m_initial)

    # Helper to finalize both manifests on any failure
    def _serialize_clusters(cl_list):
        if not cl_list:
            return []
        return [
            {
                "cluster_id": cl.cluster_id,
                "supporting_seeds": cl.supporting_seeds,
                "seed_count": cl.seed_count,
                "top_score": cl.top_score,
                "median_score": cl.median_score,
                "representative_pose_source": cl.representative_pose_source,
            }
            for cl in cl_list
        ]

    def _finalize_both_failed(err_msg: str, h_meta=None, m_meta=None, h_cl=None, m_cl=None):
        build_and_save_manifest(
            run_dir=h_run_dir,
            run_id=h_run_id,
            status="run failed",
            target_meta=dict(healthy_target_meta, pair_validation=pair_validation_meta or {}),
            binding_site_meta=binding_site_meta_healthy,
            ligand_meta=ligand_meta_initial,
            prep_meta={},
            engine_meta=engine_meta,
            grid_meta=grid_meta_healthy,
            results_meta={"clusters": _serialize_clusters(h_cl), "per_seed": h_meta or [], "failure_reason": err_msg},
            warnings=[err_msg],
            input_fingerprint=input_fingerprint,
        )
        build_and_save_manifest(
            run_dir=m_run_dir,
            run_id=m_run_id,
            status="run failed",
            target_meta=dict(mutant_target_meta, pair_validation=pair_validation_meta or {}),
            binding_site_meta=binding_site_meta_mutant,
            ligand_meta=ligand_meta_initial,
            prep_meta={},
            engine_meta=engine_meta,
            grid_meta=grid_meta_mutant,
            results_meta={"clusters": _serialize_clusters(m_cl), "per_seed": m_meta or [], "failure_reason": err_msg},
            warnings=[err_msg],
            input_fingerprint=input_fingerprint,
        )

    # 4. Receptor preparation with identical policy
    try:
        h_clean_pdb, h_rec_pdbqt, h_prep_rep = prepare_receptor(healthy_pdb_text, h_run_dir, keep_waters=keep_waters)
        m_clean_pdb, m_rec_pdbqt, m_prep_rep = prepare_receptor(mutant_pdb_text, m_run_dir, keep_waters=keep_waters)
    except Exception as e:
        _finalize_both_failed(f"Receptor preparation failure: {e}")
        raise

    # 5. Prepare ligand identically
    try:
        h_sdf, h_lig_pdbqt, h_lig_rep = prepare_ligand_from_smiles(
            ligand_smiles, ligand_name, h_run_dir, ph=7.4, pubchem_cid=pubchem_cid
        )
        m_sdf, m_lig_pdbqt, m_lig_rep = prepare_ligand_from_smiles(
            ligand_smiles, ligand_name, m_run_dir, ph=7.4, pubchem_cid=pubchem_cid
        )
    except Exception as e:
        _finalize_both_failed(f"Ligand preparation failure: {e}")
        raise

    h_lig_hash = compute_file_sha256(h_lig_pdbqt)
    m_lig_hash = compute_file_sha256(m_lig_pdbqt)

    # 6. Verify matched conditions.  Domains come from each target's own
    # metadata (never from one shared string), so a real domain mismatch is
    # detected instead of being compared with itself.
    effective_healthy_domain = healthy_domain or healthy_target_meta.get("target_domain") or target_domain
    effective_mutant_domain = mutant_domain or mutant_target_meta.get("target_domain") or target_domain
    valid, reason = validate_matched_conditions(
        healthy_domain=effective_healthy_domain,
        mutant_domain=effective_mutant_domain,
        healthy_ligand_hash=h_lig_hash,
        mutant_ligand_hash=m_lig_hash,
        healthy_waters_policy=h_prep_rep["waters_policy"],
        mutant_waters_policy=m_prep_rep["waters_policy"],
        healthy_seeds=seeds,
        mutant_seeds=seeds,
        healthy_exhaustiveness=exhaustiveness,
        mutant_exhaustiveness=exhaustiveness,
    )
    if not valid:
        err_msg = f"Matched comparison rejected: {reason}"
        _finalize_both_failed(err_msg)
        raise ValueError(err_msg)

    # 7. Execute multi-seed Vina for Healthy
    try:
        h_poses, h_meta = run_vina_multi_seeds(
            receptor_pdbqt=h_rec_pdbqt,
            ligand_pdbqt=h_lig_pdbqt,
            grid=h_grid,
            output_dir=h_run_dir,
            seeds=seeds,
            exhaustiveness=exhaustiveness,
        )
        h_clusters = cluster_poses_across_seeds(h_poses, rmsd_threshold=2.0)
    except Exception as e:
        _finalize_both_failed(f"Healthy docking engine failure: {e}")
        raise

    # 8. Execute multi-seed Vina for Mutant
    try:
        m_poses, m_meta = run_vina_multi_seeds(
            receptor_pdbqt=m_rec_pdbqt,
            ligand_pdbqt=m_lig_pdbqt,
            grid=m_grid,
            output_dir=m_run_dir,
            seeds=seeds,
            exhaustiveness=exhaustiveness,
        )
        m_clusters = cluster_poses_across_seeds(m_poses, rmsd_threshold=2.0)
    except Exception as e:
        _finalize_both_failed(f"Mutant docking engine failure: {e}", h_meta=h_meta, h_cl=h_clusters)
        raise

    # 9. Consensus cluster selection, multi-seed support and binding-mode gate
    binding_mode_details: dict = {}
    try:
        h_top_cluster, m_top_cluster, delta_score, reprod_warnings = verify_comparison_clusters(
            h_clusters,
            m_clusters,
            healthy_receptor_pdb=h_clean_pdb,
            mutant_receptor_pdb=m_clean_pdb,
            healthy_pocket_residues=healthy_pocket_residues,
            mutant_pocket_residues=mutant_pocket_residues,
            details=binding_mode_details,
        )
    except ValueError as val_err:
        _finalize_both_failed(str(val_err), h_meta=h_meta, m_meta=m_meta, h_cl=h_clusters, m_cl=m_clusters)
        raise

    # 10. Scientific validation status via optional WT reference redocking with strict chemical identity gate
    validation_status = "exploratory—no reference validation"
    validation_reason = "No verified WT co-crystal ligand was available for redocking validation."
    redocking_meta = {}
    identity_gate_meta = {}
    if wt_reference_ligand_meta and wt_reference_ligand_meta.get("ref_ligand_pdb_text"):
        comp_id = wt_reference_ligand_meta.get("resname", "")
        ref_smiles = wt_reference_ligand_meta.get("ref_smiles")
        if not ref_smiles and comp_id:
            ref_smiles = fetch_chem_comp_smiles(comp_id)

        ref_smiles_hash = hashlib.sha256(ref_smiles.encode("utf-8")).hexdigest() if ref_smiles else None
        sub_smiles_hash = hashlib.sha256(ligand_smiles.encode("utf-8")).hexdigest() if ligand_smiles else None

        id_passed, id_reason, id_details = verify_reference_ligand_identity(ref_smiles or "", ligand_smiles)
        identity_gate_meta = {
            "reference_component_id": comp_id,
            "reference_smiles_sha256": ref_smiles_hash,
            "submitted_smiles_sha256": sub_smiles_hash,
            "identity_gate_passed": id_passed,
            "reason": id_reason,
            "details": id_details,
        }

        if not id_passed:
            validation_status = "not validated for affinity comparison"
            validation_reason = (
                "The submitted ligand is not chemically identical to the verified "
                "WT co-crystal reference ligand."
            )
            reprod_warnings.append(
                f"Validation note: {validation_reason}"
            )
        else:
            ref_text = wt_reference_ligand_meta["ref_ligand_pdb_text"]
            try:
                val_details = validate_redocking_top_n(
                    ref_text,
                    h_poses,
                    h_clusters,
                    top_n=10,
                    rmsd_cutoff=2.0,
                    min_seed_support=2,
                )
                redocking_meta = val_details
                if val_details.get("protocol_validated"):
                    validation_status = "protocol validated on WT reference"
                    validation_reason = (
                        "The submitted ligand matched the verified WT reference ligand and the "
                        "top-ranked pose reproduced the crystal pose within the RMSD cutoff with "
                        "support from at least two independent seeds."
                    )
                elif val_details.get("top1_rmsd_pass") and not val_details.get("top1_cluster_pass"):
                    validation_status = "not validated for affinity comparison"
                    validation_reason = (
                        "The top-ranked pose reproduced the crystal pose, but its cluster lacked "
                        "support from at least two independent seeds."
                    )
                elif val_details.get("recovery_rmsd_pass"):
                    validation_status = "not validated for affinity comparison"
                    validation_reason = (
                        "A pose within the Top-N reproduced the crystal pose, but the top-ranked "
                        "pose did not. Top-N recovery is reported separately and is not protocol validation."
                    )
                else:
                    validation_status = "not validated for affinity comparison"
                    validation_reason = (
                        "Redocking validation did not meet the top-1 RMSD or seed-support criteria."
                    )
            except Exception as val_ex:
                validation_status = "not validated for affinity comparison"
                validation_reason = f"WT reference validation failed: {val_ex}"
                reprod_warnings.append(f"WT reference validation failed: {val_ex}")

    # A WT-reference redocking result validates the healthy arm only.  Copying
    # that status to the mutant manifest would falsely claim the mutant was
    # validated, so the mutant arm gets its own explicit status.
    mutant_validation_status = validation_status
    mutant_validation_reason = validation_reason
    if validation_status == "protocol validated on WT reference":
        mutant_validation_status = "matched comparison—mutant arm (WT validation not transferable)"
        mutant_validation_reason = (
            "The WT-reference redocking gate validates the healthy arm only. "
            "The mutant pose was not independently validated against a mutant co-crystal."
        )

    # 11. Save successful manifests for both runs
    build_and_save_manifest(
        run_dir=h_run_dir,
        run_id=h_run_id,
        status=validation_status,
        target_meta=dict(healthy_target_meta, pair_validation=pair_validation_meta or {}),
        binding_site_meta=binding_site_meta_healthy,
        ligand_meta=dict(h_lig_rep, protonation_note=PROTONATION_ASSUMPTION_NOTE),
        prep_meta=h_prep_rep,
        engine_meta=engine_meta,
        grid_meta=grid_meta_healthy,
        results_meta={
            "top_score": min(p.score for p in h_poses),
            "dominant_cluster_id": h_top_cluster.cluster_id,
            "dominant_cluster_median": h_top_cluster.median_score,
            "dominant_cluster_seed_count": h_top_cluster.seed_count,
            "clusters_count": len(h_clusters),
            "clusters": _serialize_clusters(h_clusters),
            "per_seed": h_meta,
            "validation": {
                "status": validation_status,
                "reason": validation_reason,
            },
            "binding_mode_check": binding_mode_details,
            "redocking_validation": redocking_meta,
            "reference_ligand_identity_gate": identity_gate_meta,
        },
        warnings=reprod_warnings,
        input_fingerprint=input_fingerprint,
    )
    build_and_save_manifest(
        run_dir=m_run_dir,
        run_id=m_run_id,
        status=mutant_validation_status,
        target_meta=dict(mutant_target_meta, pair_validation=pair_validation_meta or {}),
        binding_site_meta=binding_site_meta_mutant,
        ligand_meta=dict(m_lig_rep, protonation_note=PROTONATION_ASSUMPTION_NOTE),
        prep_meta=m_prep_rep,
        engine_meta=engine_meta,
        grid_meta=grid_meta_mutant,
        results_meta={
            "top_score": min(p.score for p in m_poses),
            "dominant_cluster_id": m_top_cluster.cluster_id,
            "dominant_cluster_median": m_top_cluster.median_score,
            "dominant_cluster_seed_count": m_top_cluster.seed_count,
            "clusters_count": len(m_clusters),
            "clusters": _serialize_clusters(m_clusters),
            "per_seed": m_meta,
            "validation": {
                "status": mutant_validation_status,
                "reason": mutant_validation_reason,
            },
            "binding_mode_check": binding_mode_details,
            "reference_ligand_identity_gate": identity_gate_meta,
        },
        warnings=reprod_warnings,
        input_fingerprint=input_fingerprint,
    )

    return {
        "comparison_label": COMPARISON_LABEL,
        "delta_score": delta_score,
        "validation_status": validation_status,
        "validation_reason": validation_reason,
        "mutant_validation_status": mutant_validation_status,
        "mutant_validation_reason": mutant_validation_reason,
        "binding_mode_check": binding_mode_details,
        "healthy": {
            "run_id": h_run_id,
            "target_id": healthy_target_meta.get("pdb_id"),
            "dominant_cluster_id": h_top_cluster.cluster_id,
            "dominant_cluster_top_score": h_top_cluster.top_score,
            "dominant_cluster_median": h_top_cluster.median_score,
            "dominant_cluster_seeds": h_top_cluster.supporting_seeds,
            "dominant_cluster_pose": h_top_cluster.representative_pose,
            "dominant_cluster_pose_source": h_top_cluster.representative_pose_source,
            "seed_count": h_top_cluster.seed_count,
            "clean_pdb": h_clean_pdb,
            "poses": h_poses,
            "clusters": h_clusters,
            "grid": h_grid,
            "pocket_residues": healthy_pocket_residues,
        },
        "mutant": {
            "run_id": m_run_id,
            "target_id": mutant_target_meta.get("pdb_id"),
            "dominant_cluster_id": m_top_cluster.cluster_id,
            "dominant_cluster_top_score": m_top_cluster.top_score,
            "dominant_cluster_median": m_top_cluster.median_score,
            "dominant_cluster_seeds": m_top_cluster.supporting_seeds,
            "dominant_cluster_pose": m_top_cluster.representative_pose,
            "dominant_cluster_pose_source": m_top_cluster.representative_pose_source,
            "seed_count": m_top_cluster.seed_count,
            "clean_pdb": m_clean_pdb,
            "poses": m_poses,
            "clusters": m_clusters,
            "grid": m_grid,
            "pocket_residues": mutant_pocket_residues,
        },
        "pair_validation": pair_validation_meta or {},
        "redocking_validation": redocking_meta,
        "reference_ligand_identity_gate": identity_gate_meta,
        "warnings": reprod_warnings,
        "disclaimer": DISCLAIMER_TEXT,
        "protonation_note": PROTONATION_ASSUMPTION_NOTE,
    }


def run_exploratory_cross_docking(
    receptor1_pdb_text: str,
    receptor2_pdb_text: str,
    receptor1_target_meta: dict,
    receptor2_target_meta: dict,
    ligand_smiles: str,
    ligand_name: str,
    receptor1_pocket_residues: list[str],
    receptor2_pocket_residues: list[str],
    padding: float = 8.0,
    seeds: list[int] | None = None,
    exhaustiveness: int = 8,
    input_fingerprint: str | None = None,
    pubchem_cid: int | None = None,
    keep_waters: bool = False,
) -> dict:
    """
    Execute exploratory cross-protein docking between any two structures and independently specified pockets.
    Saves complete provenance and individual manifests with status="exploratory—cross-protein docking".
    Never produces a matched delta score or comparative ranking, preserving scientific integrity.
    """
    if seeds is None:
        seeds = [42, 101, 2024]
    if not receptor1_pdb_text or not receptor2_pdb_text:
        raise ValueError("Missing coordinates for Receptor 1 or Receptor 2.")
    if not receptor1_pocket_residues or not receptor2_pocket_residues:
        raise ValueError("Exploratory docking requires pocket residues defined for both receptors.")

    if receptor1_target_meta.get("selected_chain"):
        validate_pocket_chain(receptor1_pocket_residues, receptor1_target_meta["selected_chain"], "Receptor 1")
    if receptor2_target_meta.get("selected_chain"):
        validate_pocket_chain(receptor2_pocket_residues, receptor2_target_meta["selected_chain"], "Receptor 2")

    vina_ver = get_vina_version()
    r1_created_at = receptor1_target_meta.get("created_at_utc") or datetime.datetime.now(datetime.timezone.utc).isoformat()
    r2_created_at = receptor2_target_meta.get("created_at_utc") or datetime.datetime.now(datetime.timezone.utc).isoformat()

    # 1. Setup independent run directories
    r1_run_id, r1_run_dir = create_run_directory("dock_cross_receptor1")
    r2_run_id, r2_run_dir = create_run_directory("dock_cross_receptor2")

    # Save receptor_source.pdb immediately
    with open(r1_run_dir / "receptor_source.pdb", "w", encoding="utf-8") as f:
        f.write(receptor1_pdb_text)
    with open(r2_run_dir / "receptor_source.pdb", "w", encoding="utf-8") as f:
        f.write(receptor2_pdb_text)

    # 2. Compute individual pocket grids
    r1_grid = define_grid_from_residues(receptor1_pdb_text, receptor1_pocket_residues, padding=padding)
    r2_grid = define_grid_from_residues(receptor2_pdb_text, receptor2_pocket_residues, padding=padding)

    grid_meta_r1 = {
        "center": [r1_grid.center_x, r1_grid.center_y, r1_grid.center_z],
        "size": [r1_grid.size_x, r1_grid.size_y, r1_grid.size_z],
    }
    grid_meta_r2 = {
        "center": [r2_grid.center_x, r2_grid.center_y, r2_grid.center_z],
        "size": [r2_grid.size_x, r2_grid.size_y, r2_grid.size_z],
    }

    engine_meta = {
        "name": "AutoDock Vina",
        "version": vina_ver,
        "seeds": seeds,
        "exhaustiveness": exhaustiveness,
        "num_modes": 9,
        "energy_range": 3.0,
    }

    binding_site_meta_r1 = {
        "method": "residues",
        "residues": receptor1_pocket_residues,
        "padding": padding,
    }
    binding_site_meta_r2 = {
        "method": "residues",
        "residues": receptor2_pocket_residues,
        "padding": padding,
    }

    ligand_meta_initial = {
        "name": ligand_name,
        "input_smiles": ligand_smiles,
        "pubchem_cid": pubchem_cid,
        "protonation_note": PROTONATION_ASSUMPTION_NOTE,
    }

    # 3. Create initial manifest for both runs BEFORE preparation or engine runs
    r1_initial = {
        "schema_version": 1,
        "run_id": r1_run_id,
        "status": "running",
        "created_at_utc": r1_created_at,
        "target": receptor1_target_meta,
        "binding_site": binding_site_meta_r1,
        "ligand": ligand_meta_initial,
        "grid": grid_meta_r1,
        "engine": engine_meta,
    }
    write_initial_manifest(r1_run_dir, r1_initial)

    r2_initial = {
        "schema_version": 1,
        "run_id": r2_run_id,
        "status": "running",
        "created_at_utc": r2_created_at,
        "target": receptor2_target_meta,
        "binding_site": binding_site_meta_r2,
        "ligand": ligand_meta_initial,
        "grid": grid_meta_r2,
        "engine": engine_meta,
    }
    write_initial_manifest(r2_run_dir, r2_initial)

    def _serialize_clusters(cl_list):
        if not cl_list:
            return []
        return [
            {
                "cluster_id": cl.cluster_id,
                "supporting_seeds": cl.supporting_seeds,
                "seed_count": cl.seed_count,
                "top_score": cl.top_score,
                "median_score": cl.median_score,
                "representative_pose_source": cl.representative_pose_source,
            }
            for cl in cl_list
        ]

    def _finalize_cross_failed(err_msg: str, r1_m=None, r2_m=None, r1_cl=None, r2_cl=None, r1_rep=None, r2_rep=None, r1_lrep=None, r2_lrep=None):
        build_and_save_manifest(
            run_dir=r1_run_dir,
            run_id=r1_run_id,
            status="run failed",
            target_meta=receptor1_target_meta,
            binding_site_meta=binding_site_meta_r1,
            ligand_meta=dict(r1_lrep or ligand_meta_initial, protonation_note=PROTONATION_ASSUMPTION_NOTE),
            prep_meta=r1_rep or {},
            engine_meta=engine_meta,
            grid_meta=grid_meta_r1,
            results_meta={"clusters": _serialize_clusters(r1_cl), "per_seed": r1_m or [], "failure_reason": err_msg},
            warnings=[err_msg],
            input_fingerprint=input_fingerprint,
        )
        build_and_save_manifest(
            run_dir=r2_run_dir,
            run_id=r2_run_id,
            status="run failed",
            target_meta=receptor2_target_meta,
            binding_site_meta=binding_site_meta_r2,
            ligand_meta=dict(r2_lrep or ligand_meta_initial, protonation_note=PROTONATION_ASSUMPTION_NOTE),
            prep_meta=r2_rep or {},
            engine_meta=engine_meta,
            grid_meta=grid_meta_r2,
            results_meta={"clusters": _serialize_clusters(r2_cl), "per_seed": r2_m or [], "failure_reason": err_msg},
            warnings=[err_msg],
            input_fingerprint=input_fingerprint,
        )

    # 4. Preparation
    try:
        r1_clean_pdb, r1_rec_pdbqt, r1_prep_rep = prepare_receptor(receptor1_pdb_text, r1_run_dir, keep_waters=keep_waters)
        r2_clean_pdb, r2_rec_pdbqt, r2_prep_rep = prepare_receptor(receptor2_pdb_text, r2_run_dir, keep_waters=keep_waters)
    except Exception as e:
        _finalize_cross_failed(f"Receptor preparation failure: {e}")
        raise

    try:
        r1_sdf, r1_lig_pdbqt, r1_lig_rep = prepare_ligand_from_smiles(
            ligand_smiles, ligand_name, r1_run_dir, ph=7.4, pubchem_cid=pubchem_cid
        )
        r2_sdf, r2_lig_pdbqt, r2_lig_rep = prepare_ligand_from_smiles(
            ligand_smiles, ligand_name, r2_run_dir, ph=7.4, pubchem_cid=pubchem_cid
        )
    except Exception as e:
        _finalize_cross_failed(f"Ligand preparation failure: {e}", r1_rep=r1_prep_rep, r2_rep=r2_prep_rep)
        raise

    # 5. Run Vina for Receptor 1
    try:
        r1_poses, r1_meta = run_vina_multi_seeds(
            receptor_pdbqt=r1_rec_pdbqt,
            ligand_pdbqt=r1_lig_pdbqt,
            grid=r1_grid,
            output_dir=r1_run_dir,
            seeds=seeds,
            exhaustiveness=exhaustiveness,
        )
        r1_clusters = cluster_poses_across_seeds(r1_poses, rmsd_threshold=2.0)
    except Exception as e:
        _finalize_cross_failed(f"Receptor 1 docking failure: {e}", r1_rep=r1_prep_rep, r2_rep=r2_prep_rep, r1_lrep=r1_lig_rep, r2_lrep=r2_lig_rep)
        raise

    # 6. Run Vina for Receptor 2
    try:
        r2_poses, r2_meta = run_vina_multi_seeds(
            receptor_pdbqt=r2_rec_pdbqt,
            ligand_pdbqt=r2_lig_pdbqt,
            grid=r2_grid,
            output_dir=r2_run_dir,
            seeds=seeds,
            exhaustiveness=exhaustiveness,
        )
        r2_clusters = cluster_poses_across_seeds(r2_poses, rmsd_threshold=2.0)
    except Exception as e:
        _finalize_cross_failed(f"Receptor 2 docking failure: {e}", r1_m=r1_meta, r1_cl=r1_clusters, r1_rep=r1_prep_rep, r2_rep=r2_prep_rep, r1_lrep=r1_lig_rep, r2_lrep=r2_lig_rep)
        raise

    # 7. Consensus cluster selection and reproducibility checks
    try:
        r1_top_cluster = _select_consensus_cluster(r1_clusters)
        r2_top_cluster = _select_consensus_cluster(r2_clusters)
    except Exception as e:
        _finalize_cross_failed(f"Cluster selection failure: {e}", r1_m=r1_meta, r2_m=r2_meta, r1_cl=r1_clusters, r2_cl=r2_clusters, r1_rep=r1_prep_rep, r2_rep=r2_prep_rep, r1_lrep=r1_lig_rep, r2_lrep=r2_lig_rep)
        raise

    r1_cluster_reproducible = r1_top_cluster.seed_count >= MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT
    r2_cluster_reproducible = r2_top_cluster.seed_count >= MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT

    # 8. Build and save manifests with status="exploratory—cross-protein docking"
    CROSS_STATUS = "exploratory—cross-protein docking"
    build_and_save_manifest(
        run_dir=r1_run_dir,
        run_id=r1_run_id,
        status=CROSS_STATUS,
        target_meta=receptor1_target_meta,
        binding_site_meta=binding_site_meta_r1,
        ligand_meta=dict(r1_lig_rep, protonation_note=PROTONATION_ASSUMPTION_NOTE),
        prep_meta=r1_prep_rep,
        engine_meta=engine_meta,
        grid_meta=grid_meta_r1,
        results_meta={
            "top_score": min(p.score for p in r1_poses),
            "dominant_cluster_id": r1_top_cluster.cluster_id,
            "dominant_cluster_median": r1_top_cluster.median_score,
            "dominant_cluster_seed_count": r1_top_cluster.seed_count,
            "cluster_reproducible": r1_cluster_reproducible,
            "clusters": _serialize_clusters(r1_clusters),
            "per_seed": r1_meta,
        },
        warnings=[],
        input_fingerprint=input_fingerprint,
    )

    build_and_save_manifest(
        run_dir=r2_run_dir,
        run_id=r2_run_id,
        status=CROSS_STATUS,
        target_meta=receptor2_target_meta,
        binding_site_meta=binding_site_meta_r2,
        ligand_meta=dict(r2_lig_rep, protonation_note=PROTONATION_ASSUMPTION_NOTE),
        prep_meta=r2_prep_rep,
        engine_meta=engine_meta,
        grid_meta=grid_meta_r2,
        results_meta={
            "top_score": min(p.score for p in r2_poses),
            "dominant_cluster_id": r2_top_cluster.cluster_id,
            "dominant_cluster_median": r2_top_cluster.median_score,
            "dominant_cluster_seed_count": r2_top_cluster.seed_count,
            "cluster_reproducible": r2_cluster_reproducible,
            "clusters": _serialize_clusters(r2_clusters),
            "per_seed": r2_meta,
        },
        warnings=[],
        input_fingerprint=input_fingerprint,
    )

    EXPLORATORY_CROSS_LABEL = (
        "Exploratory Cross-Protein Docking: Independent separate simulations. "
        "No matched delta score or comparative mutation claim is permitted between different proteins."
    )

    return {
        "status": CROSS_STATUS,
        "comparison_label": EXPLORATORY_CROSS_LABEL,
        "validation_status": CROSS_STATUS,
        "receptor1": {
            "run_id": r1_run_id,
            "target_id": receptor1_target_meta.get("pdb_id"),
            "dominant_cluster_id": r1_top_cluster.cluster_id,
            "dominant_cluster_top_score": r1_top_cluster.top_score,
            "dominant_cluster_median": r1_top_cluster.median_score,
            "dominant_cluster_seeds": r1_top_cluster.supporting_seeds,
            "dominant_cluster_pose": r1_top_cluster.representative_pose,
            "dominant_cluster_pose_source": r1_top_cluster.representative_pose_source,
            "seed_count": r1_top_cluster.seed_count,
            "cluster_reproducible": r1_cluster_reproducible,
            "dominant_cluster_seed_count": r1_top_cluster.seed_count,
            "clean_pdb": r1_clean_pdb,
            "poses": r1_poses,
            "clusters": r1_clusters,
            "grid": r1_grid,
            "pocket_residues": receptor1_pocket_residues,
        },
        "receptor2": {
            "run_id": r2_run_id,
            "target_id": receptor2_target_meta.get("pdb_id"),
            "dominant_cluster_id": r2_top_cluster.cluster_id,
            "dominant_cluster_top_score": r2_top_cluster.top_score,
            "dominant_cluster_median": r2_top_cluster.median_score,
            "dominant_cluster_seeds": r2_top_cluster.supporting_seeds,
            "dominant_cluster_pose": r2_top_cluster.representative_pose,
            "dominant_cluster_pose_source": r2_top_cluster.representative_pose_source,
            "seed_count": r2_top_cluster.seed_count,
            "cluster_reproducible": r2_cluster_reproducible,
            "dominant_cluster_seed_count": r2_top_cluster.seed_count,
            "clean_pdb": r2_clean_pdb,
            "poses": r2_poses,
            "clusters": r2_clusters,
            "grid": r2_grid,
            "pocket_residues": receptor2_pocket_residues,
        },
        "disclaimer": DISCLAIMER_TEXT,
        "protonation_note": PROTONATION_ASSUMPTION_NOTE,
    }
