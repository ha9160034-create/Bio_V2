# Docking validation, atom mapping, and RMSD calculations
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from docking.openbabel_io import read_text, write_text
from docking.models import DockingPose, PoseCluster


def verify_reference_ligand_identity(
    reference_smiles: str,
    submitted_smiles: str,
) -> tuple[bool, str, dict]:
    """Strictly verify that the submitted ligand is chemically identical to the WT reference ligand.

    Standardizes both molecules and checks:
    1. canonical isomeric SMILES match
    2. heavy atom count match
    3. Hill system molecular formula match
    """
    if not (reference_smiles or "").strip() or not (submitted_smiles or "").strip():
        return False, "Missing reference or submitted SMILES for identity gate.", {}

    mol_ref = Chem.MolFromSmiles(reference_smiles.strip())
    mol_sub = Chem.MolFromSmiles(submitted_smiles.strip())

    if mol_ref is None or mol_sub is None:
        return False, "Failed to parse reference or submitted SMILES with RDKit.", {}

    # Documented standardization: canonical isomeric SMILES with explicit stereo and aromaticity
    ref_can = Chem.MolToSmiles(mol_ref, isomericSmiles=True, canonical=True)
    sub_can = Chem.MolToSmiles(mol_sub, isomericSmiles=True, canonical=True)

    ref_heavy = mol_ref.GetNumHeavyAtoms()
    sub_heavy = mol_sub.GetNumHeavyAtoms()

    ref_formula = CalcMolFormula(mol_ref)
    sub_formula = CalcMolFormula(mol_sub)

    details = {
        "reference_canonical_smiles": ref_can,
        "submitted_canonical_smiles": sub_can,
        "reference_heavy_atom_count": ref_heavy,
        "submitted_heavy_atom_count": sub_heavy,
        "reference_molecular_formula": ref_formula,
        "submitted_molecular_formula": sub_formula,
    }

    if ref_can != sub_can:
        return False, f"Canonical isomeric SMILES mismatch: '{ref_can}' != '{sub_can}'.", details
    if ref_heavy != sub_heavy:
        return False, f"Heavy atom count mismatch: {ref_heavy} != {sub_heavy}.", details
    if ref_formula != sub_formula:
        return False, f"Molecular formula mismatch: '{ref_formula}' != '{sub_formula}'.", details

    return True, "Chemical identity verified.", details


def parse_pdbqt_coords(pdbqt_block: str) -> np.ndarray:
    """Extract coordinates of heavy atoms from a PDBQT block (first model)."""
    coords = []
    for line in pdbqt_block.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_type = line[77:].strip() if len(line) >= 78 else ""
            atom_name = line[12:16].strip() if len(line) >= 16 else ""
            is_h = atom_type.upper().startswith("H") or (not atom_type and atom_name.upper().startswith("H"))
            if not is_h:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass
    return np.array(coords, dtype=np.float32)

def calculate_heavy_atom_rmsd(ref_coords: np.ndarray, pred_coords: np.ndarray) -> float:
    """Calculate heavy-atom RMSD requiring strictly identical atom counts."""
    if len(ref_coords) != len(pred_coords) or len(ref_coords) == 0:
        raise ValueError(
            f"Cannot calculate RMSD: unequal atom count (ref={len(ref_coords)}, pred={len(pred_coords)})."
        )
    diff = ref_coords - pred_coords
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))


def pose_block_to_mol(pdbqt_block: str) -> Chem.Mol | None:
    """Convert a PDBQT pose block into an RDKit molecule (heavy atoms only).

    PDBQT carries no bond orders, so OpenBabel perceives the chemical graph and
    the result is kept unsanitised.  Returns ``None`` when conversion fails.
    """
    try:
        ob_pose = read_text(pdbqt_block, "pdbqt")
        mol = Chem.MolFromPDBBlock(write_text(ob_pose, "pdb"), sanitize=False)
        if mol is None:
            mol = Chem.MolFromMolBlock(write_text(ob_pose, "mol"), sanitize=False)
        if mol is None:
            return None
        mol = Chem.RemoveHs(mol, sanitize=False)
        if mol.GetNumConformers() == 0:
            return None
        return mol
    except Exception:
        return None


def symmetric_pose_rmsd(block_a: str, block_b: str) -> tuple[float, str]:
    """No-fit, symmetry-aware heavy-atom RMSD between two docking poses.

    Uses RDKit's ``CalcRMS`` (the documented function for comparing docking
    poses without superposition) so that chemically equivalent atoms in
    symmetric groups are not counted as different binding modes.  Falls back to
    a raw coordinate RMSD only when the chemical graph cannot be perceived,
    and reports which method was used.
    """
    raw_a = parse_pdbqt_coords(block_a)
    raw_b = parse_pdbqt_coords(block_b)
    if len(raw_a) != len(raw_b) or len(raw_a) == 0:
        raise ValueError(
            f"Cannot compare poses: unequal heavy-atom counts ({len(raw_a)} vs {len(raw_b)})."
        )
    raw_rmsd = float(np.sqrt(np.mean(np.sum((raw_a - raw_b) ** 2, axis=1))))

    mol_a = pose_block_to_mol(block_a)
    mol_b = pose_block_to_mol(block_b)
    if mol_a is None or mol_b is None or mol_a.GetNumAtoms() != mol_b.GetNumAtoms():
        return round(raw_rmsd, 3), "heavy_atom_euclidean_fallback"
    try:
        rmsd = float(rdMolAlign.CalcRMS(mol_a, mol_b))
        return round(rmsd, 3), "rdkit_symmetry_aware_calc_rms"
    except Exception:
        return round(raw_rmsd, 3), "heavy_atom_euclidean_fallback"

def calculate_mapped_redocking_rmsd(ref_ligand_pdb_text: str, pose_pdbqt_block: str) -> tuple[float, dict]:
    """
    Calculate heavy-atom RMSD between reference crystal ligand and a docked pose
    using Maximum Common Substructure (MCS) topology-invariant atom mapping.
    Evaluates every pair of reference and pose mappings to handle chemical symmetry correctly.
    Refuses incomplete mappings or unequal heavy-atom counts.
    Do not fit or rotate the ligand before RMSD; the docked and crystal ligand share the receptor coordinate frame.
    """
    if not ref_ligand_pdb_text.strip() or not pose_pdbqt_block.strip():
        raise ValueError("Missing reference ligand or pose coordinates.")

    # 1. Parse reference ligand into RDKit molecule
    mol_ref = Chem.MolFromPDBBlock(ref_ligand_pdb_text, sanitize=False)
    if mol_ref is None:
        ob_ref = read_text(ref_ligand_pdb_text, "pdb")
        mol_ref = Chem.MolFromMolBlock(write_text(ob_ref, "mol"), sanitize=False)

    if mol_ref is None:
        raise ValueError("Could not parse reference ligand coordinates into a chemical molecule.")

    # 2. Convert pose PDBQT block into RDKit molecule via OpenBabel to preserve coordinates
    ob_pose = read_text(pose_pdbqt_block, "pdbqt")
    mol_pose = Chem.MolFromPDBBlock(write_text(ob_pose, "pdb"), sanitize=False)
    if mol_pose is None:
        mol_pose = Chem.MolFromMolBlock(write_text(ob_pose, "mol"), sanitize=False)

    if mol_pose is None:
        raise ValueError("Could not parse docked pose coordinates into a chemical molecule.")

    # Strip explicit/polar hydrogens to strictly evaluate heavy-atom RMSD while preserving conformer coordinates
    mol_ref = Chem.RemoveHs(mol_ref, sanitize=False)
    mol_pose = Chem.RemoveHs(mol_pose, sanitize=False)

    ref_heavy_count = mol_ref.GetNumHeavyAtoms()
    pose_heavy_count = mol_pose.GetNumHeavyAtoms()

    if ref_heavy_count != pose_heavy_count or ref_heavy_count == 0:
        raise ValueError(
            f"Validation unavailable: reference ligand has {ref_heavy_count} heavy atoms, "
            f"docked ligand has {pose_heavy_count} heavy atoms. Molecules must match exactly."
        )

    if mol_ref.GetNumConformers() == 0 or mol_pose.GetNumConformers() == 0:
        raise ValueError("Validation unavailable: ligand molecule does not contain 3D conformer coordinates.")

    Chem.FastFindRings(mol_ref)
    Chem.FastFindRings(mol_pose)

    # 3. Establish documented atom mapping via MCS.
    # Atom mapping always uses CompareAny because PDBQT carries no bond orders
    # and OpenBabel may assign a different Kekule form for the same aromatic
    # ring; restricting the mapping to CompareOrder would hide valid symmetry
    # matches.  A separate CompareOrder MCS is computed as a diagnostic of
    # whether a fully bond-order-consistent substructure exists.
    common_kwargs = dict(
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        timeout=10,
    )
    mcs = rdFMCS.FindMCS(
        [mol_ref, mol_pose],
        bondCompare=rdFMCS.BondCompare.CompareAny,
        **common_kwargs,
    )
    strict_mcs = rdFMCS.FindMCS(
        [mol_ref, mol_pose],
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        **common_kwargs,
    )
    bond_order_strict = strict_mcs.numAtoms == ref_heavy_count and ref_heavy_count > 0
    mapping_method = "rdkit_mcs_bond_order_agnostic_exhaustive_pairs"

    if mcs.numAtoms != ref_heavy_count:
        raise ValueError(
            f"Validation unavailable: atom mapping could only match {mcs.numAtoms}/{ref_heavy_count} heavy atoms. "
            "Chemical identity mismatch between reference and docked ligand."
        )

    query_mol = Chem.MolFromSmarts(mcs.smartsString)
    match_ref_all = mol_ref.GetSubstructMatches(query_mol, uniquify=False)
    match_pose_all = mol_pose.GetSubstructMatches(query_mol, uniquify=False)

    if not match_ref_all or not match_pose_all:
        raise ValueError("Validation unavailable: could not map atoms from MCS query.")

    conf_ref = mol_ref.GetConformer()
    conf_pose = mol_pose.GetConformer()

    # 4. Iterate through every pair of reference and pose mappings to find minimum RMSD
    best_rmsd = 999.0
    best_mapping = None
    mappings_evaluated = 0

    for m_ref in match_ref_all:
        for m_pose in match_pose_all:
            mappings_evaluated += 1
            diffs = []
            mapping_pairs = []
            for idx_q in range(len(m_ref)):
                idx_r = m_ref[idx_q]
                idx_p = m_pose[idx_q]
                p_r = conf_ref.GetAtomPosition(idx_r)
                p_p = conf_pose.GetAtomPosition(idx_p)
                diffs.append([p_r.x - p_p.x, p_r.y - p_p.y, p_r.z - p_p.z])
                mapping_pairs.append({"ref_atom_idx": int(idx_r), "pose_atom_idx": int(idx_p)})

            diffs_arr = np.array(diffs, dtype=np.float32)
            rmsd = float(np.sqrt(np.mean(np.sum(diffs_arr**2, axis=1))))
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_mapping = mapping_pairs

    details = {
        "mapping_method": mapping_method,
        "bond_order_strict": bond_order_strict,
        "bond_order_note": (
            "MCS required identical bond orders."
            if bond_order_strict
            else "PDBQT carries no bond orders; mapping is connectivity/element based (CompareAny)."
        ),
        "reference_heavy_atom_count": ref_heavy_count,
        "pose_heavy_atom_count": pose_heavy_count,
        "matched_atom_count": mcs.numAtoms,
        "mappings_evaluated": mappings_evaluated,
        "selected_mapping_pairs": best_mapping,
        "final_rmsd": round(best_rmsd, 3)
    }

    # Inspect the bonds implied by the selected atom mapping.  PDBQT cannot
    # store bond orders, so this is a diagnostic rather than a hard rejection;
    # it makes any ketone/enol or amide/imidic mapping visible in the manifest.
    bond_mismatches = 0
    missing_bonds = 0
    if best_mapping:
        ref_to_pose = {
            pair["ref_atom_idx"]: pair["pose_atom_idx"] for pair in best_mapping
        }
        for bond in mol_ref.GetBonds():
            a = bond.GetBeginAtomIdx()
            b = bond.GetEndAtomIdx()
            if a not in ref_to_pose or b not in ref_to_pose:
                continue
            pose_bond = mol_pose.GetBondBetweenAtoms(ref_to_pose[a], ref_to_pose[b])
            if pose_bond is None:
                missing_bonds += 1
            elif bond.GetBondType() != pose_bond.GetBondType():
                # Aromatic Kekule differences are expected after PDBQT
                # round-tripping and are not counted as errors.
                if not (bond.GetIsAromatic() and pose_bond.GetIsAromatic()):
                    bond_mismatches += 1
    details["mapped_bonds_missing_in_pose"] = missing_bonds
    details["bond_order_mismatches"] = bond_mismatches
    details["mapping_quality"] = (
        "strict_bond_orders" if bond_order_strict and bond_mismatches == 0 and missing_bonds == 0
        else "connectivity_only_review_required"
    )
    return round(best_rmsd, 3), details

def cluster_poses_across_seeds(poses: list[DockingPose], rmsd_threshold: float = 2.0) -> list[PoseCluster]:
    """
    Cluster poses across independent seeds using heavy-atom RMSD.
    Returns clustered groups with seed counts, median scores, and representative poses.
    """
    if not poses:
        return []

    # Sort poses by binding score (ascending, best/lowest first)
    sorted_poses = sorted(poses, key=lambda p: p.score)
    clusters: list[list[DockingPose]] = []
    mol_cache: dict[int, Chem.Mol | None] = {}

    def _mol(pose: DockingPose):
        key = id(pose)
        if key not in mol_cache:
            mol_cache[key] = pose_block_to_mol(pose.pdbqt_block)
        return mol_cache[key]

    for pose in sorted_poses:
        p_coords = parse_pdbqt_coords(pose.pdbqt_block)
        assigned = False
        for cl in clusters:
            rep = cl[0]
            rep_coords = parse_pdbqt_coords(rep.pdbqt_block)
            if len(p_coords) == len(rep_coords) and len(p_coords) > 0:
                # Symmetry-aware no-fit RMSD (RDKit CalcRMS) prevents equivalent
                # atoms of symmetric groups from being split into fake clusters.
                mol_p, mol_rep = _mol(pose), _mol(rep)
                if mol_p is not None and mol_rep is not None and mol_p.GetNumAtoms() == mol_rep.GetNumAtoms():
                    try:
                        rmsd = float(rdMolAlign.CalcRMS(mol_p, mol_rep))
                    except Exception:
                        rmsd = float(np.sqrt(np.mean(np.sum((p_coords - rep_coords) ** 2, axis=1))))
                else:
                    rmsd = float(np.sqrt(np.mean(np.sum((p_coords - rep_coords) ** 2, axis=1))))
                if rmsd <= rmsd_threshold:
                    cl.append(pose)
                    assigned = True
                    break
        if not assigned:
            clusters.append([pose])

    result_clusters = []
    for idx, cl in enumerate(clusters, start=1):
        supporting_seeds = sorted(list(set(p.seed for p in cl)))
        scores = [p.score for p in cl]
        median_score = float(np.median(scores))
        top_score = float(min(scores))
        for p in cl:
            p.cluster_id = idx
        result_clusters.append(PoseCluster(
            cluster_id=idx,
            seed_count=len(supporting_seeds),
            supporting_seeds=supporting_seeds,
            top_score=round(top_score, 2),
            median_score=round(median_score, 2),
            representative_pose=cl[0],
            representative_pose_source={"seed": cl[0].seed, "rank": cl[0].rank}
        ))

    return result_clusters

def validate_redocking_top_n(
    reference_ligand_pdb_text: str,
    poses: list[DockingPose],
    clusters: list[PoseCluster],
    top_n: int = 10,
    rmsd_cutoff: float = 2.0,
    min_seed_support: int = 2,
) -> dict:
    if top_n < 1:
        raise ValueError("top_n must be >= 1")
    if not poses:
        raise ValueError("No docking poses available for redocking validation.")

    ranked = sorted(poses, key=lambda p: p.score)
    top_pose = ranked[0]
    top1_rmsd, top1_mapping = calculate_mapped_redocking_rmsd(
        reference_ligand_pdb_text, top_pose.pdbqt_block
    )
    top_pose.rmsd_to_reference = top1_rmsd

    cluster_by_id = {cluster.cluster_id: cluster for cluster in clusters}
    top1_cluster = cluster_by_id.get(top_pose.cluster_id)
    top1_seed_count = top1_cluster.seed_count if top1_cluster else 0
    top1_rmsd_pass = top1_rmsd <= rmsd_cutoff
    top1_cluster_pass = top1_seed_count >= min_seed_support

    candidates = []
    for pose in ranked[:top_n]:
        rmsd, mapping = calculate_mapped_redocking_rmsd(
            reference_ligand_pdb_text, pose.pdbqt_block
        )
        pose.rmsd_to_reference = rmsd
        candidates.append((pose, rmsd, mapping))

    recovered_pose, recovered_rmsd, recovered_mapping = min(
        candidates, key=lambda item: item[1]
    )
    recovered_cluster = cluster_by_id.get(recovered_pose.cluster_id)
    seed_count = recovered_cluster.seed_count if recovered_cluster else 0

    return {
        "top_n": min(top_n, len(ranked)),
        "rmsd_cutoff": rmsd_cutoff,
        "min_seed_support": min_seed_support,
        "method": "symmetry-aware heavy-atom mapped RMSD (no-fit)",
        "top1": {
            "seed": top_pose.seed,
            "rank": top_pose.rank,
            "score": top_pose.score,
            "cluster_id": top_pose.cluster_id,
            "rmsd_angstrom": top1_rmsd,
            "mapping_details": top1_mapping,
            "cluster_seed_count": top1_seed_count,
            "supporting_seeds": top1_cluster.supporting_seeds if top1_cluster else [],
        },
        "recovery": {
            "seed": recovered_pose.seed,
            "rank": recovered_pose.rank,
            "score": recovered_pose.score,
            "cluster_id": recovered_pose.cluster_id,
            "rmsd_angstrom": recovered_rmsd,
            "cluster_seed_count": seed_count,
            "supporting_seeds": (
                recovered_cluster.supporting_seeds if recovered_cluster else []
            ),
            "mapping_details": recovered_mapping,
        },
        "candidate_count": len(candidates),
        "top1_rmsd_pass": top1_rmsd_pass,
        "top1_cluster_seed_count": top1_seed_count,
        "top1_cluster_pass": top1_cluster_pass,
        "recovery_rmsd_pass": recovered_rmsd <= rmsd_cutoff,
        "recovery_cluster_pass": seed_count >= min_seed_support,
        # A protocol is only "validated" when the top-ranked pose itself is
        # correct and reproducible.  Top-N recovery remains a separate metric.
        "protocol_validated": bool(top1_rmsd_pass and top1_cluster_pass),
        "protocol_validation_rule": (
            "top-1 pose RMSD <= cutoff AND top-1 pose cluster supported by "
            f">= {min_seed_support} independent seeds"
        ),
    }
