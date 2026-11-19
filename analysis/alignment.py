# Sequence alignment and mutation biochemical impact classification
from Bio.Align import PairwiseAligner
import streamlit as st
from analysis.constants import AA_PROPS, AA_1TO3, AA_VOLUMES

def analyze_impact(h_res, m_res, h_sasa, m_sasa):
    """Analyze the biochemical and structural impact of a mutation based on physicochemical properties and SASA (Solvent-Accessible Surface Area)."""
    # Handle insertions and deletions (indels) in a biologically accurate manner
    if h_res == '-' and m_res != '-':
        return "Residue Insertion (Critical Structural Shift)"
    if h_res != '-' and m_res == '-':
        return "Residue Deletion (Critical Structural Shift)"
    
    # Convert 1-letter codes to 3-letter codes when provided
    h_code = AA_1TO3.get(h_res.upper(), h_res.upper()) if len(h_res) == 1 else h_res.upper()
    m_code = AA_1TO3.get(m_res.upper(), m_res.upper()) if len(m_res) == 1 else m_res.upper()

    # Identical residues indicate no mutation
    if h_code == m_code:
        return "Identical (No Mutation)"
    
    h_type = AA_PROPS.get(h_code, 'Unknown')
    m_type = AA_PROPS.get(m_code, 'Unknown')
    
    impacts = []
    
    # 1. Charge and chemical class analysis
    is_h_acidic = "Acidic" in h_type
    is_m_acidic = "Acidic" in m_type
    is_h_basic = "Basic" in h_type
    is_m_basic = "Basic" in m_type
    is_h_charged = is_h_acidic or is_h_basic
    is_m_charged = is_m_acidic or is_m_basic

    if (is_h_acidic and is_m_basic) or (is_h_basic and is_m_acidic):
        impacts.append("Charge Flip (Critical)")
    elif is_h_charged and not is_m_charged:
        impacts.append("Charge Loss (Electrostatic Disruption)")
    elif not is_h_charged and is_m_charged:
        impacts.append("Charge Gain (Electrostatic Shift)")
    elif h_type != m_type:
        impacts.append(f"Chem-Class Change ({h_type} ➔ {m_type})")

    # 2. Aromaticity change (pi-stacking)
    is_h_arom = "Aromatic" in h_type
    is_m_arom = "Aromatic" in m_type
    if is_h_arom != is_m_arom:
        impacts.append("Aromaticity Change (Pi-Stacking Risk)")
    
    # 3. Disulfide bond risk
    if h_code == 'CYS' or m_code == 'CYS':
        impacts.append("Disulfide Bond Risk (Critical)")
    
    # 4. Proline substitution (helix breaker)
    if (h_code == 'PRO') != (m_code == 'PRO'):
        impacts.append("Proline Substitution (Helix Breaker)")
    
    # 5. Glycine flexibility
    if (h_code == 'GLY') != (m_code == 'GLY'):
        impacts.append("Glycine Flexibility Shift")

    # 6. Steric volume change
    h_vol = AA_VOLUMES.get(h_code)
    m_vol = AA_VOLUMES.get(m_code)
    if h_vol is not None and m_vol is not None:
        delta_vol = m_vol - h_vol
        if delta_vol >= 55.0:
            impacts.append(f"Steric Bulk Increase (+{delta_vol:.0f} Å³ [cubic Angstroms] Clash Risk)")
        elif delta_vol <= -55.0:
            impacts.append(f"Steric Void Creation ({delta_vol:.0f} Å³ [cubic Angstroms])")

    # 7. SASA (Solvent-Accessible Surface Area) and structural environment
    try:
        diff_sasa = m_sasa - h_sasa
        if h_sasa < 15.0:
            # Buried core residue
            if is_m_charged or (h_vol and m_vol and abs(m_vol - h_vol) >= 50):
                impacts.append("Buried Core Mutation (High Destabilization Risk)")
            else:
                impacts.append("Core Mutation")

        if diff_sasa >= 15.0:
            impacts.append("Increased Exposure (+ΔSASA — Solvent-Accessible Surface Area)")
        elif diff_sasa <= -15.0:
            impacts.append("Increased Burial (-ΔSASA — Solvent-Accessible Surface Area)")
    except Exception:
        pass
        
    if not impacts:
        if h_type == m_type:
            return "Conservative Substitution (Same Chemical Class)"
        return "Minor Physicochemical Change"

    return " | ".join(impacts)

def calculate_sequence_identity(aligned_a: str, aligned_b: str) -> dict:
    """Compute sequence identity with one documented denominator.

    Identity is defined over *aligned non-gap columns only* (the same rule used
    by :mod:`docking.pair_validation` for the matched-docking >= 90 % gate).
    Gaps are reported separately so that indel-rich alignments are not silently
    penalised twice (once in the numerator and once in the denominator).

    Returns a dict with ``identity_pct``, ``identical``, ``aligned_non_gap``,
    ``gaps`` and ``alignment_length``.
    """
    if not aligned_a or not aligned_b or len(aligned_a) != len(aligned_b):
        return {
            "identity_pct": 0.0,
            "identical": 0,
            "aligned_non_gap": 0,
            "gaps": 0,
            "alignment_length": max(len(aligned_a or ""), len(aligned_b or "")),
            "valid": False,
        }

    identical = 0
    aligned_non_gap = 0
    gaps = 0
    for a, b in zip(aligned_a, aligned_b):
        if a == "-" or b == "-":
            gaps += 1
            continue
        aligned_non_gap += 1
        if a.upper() == b.upper():
            identical += 1

    identity_pct = (identical / aligned_non_gap * 100.0) if aligned_non_gap else 0.0
    return {
        "identity_pct": round(identity_pct, 2),
        "identical": identical,
        "aligned_non_gap": aligned_non_gap,
        "gaps": gaps,
        "alignment_length": len(aligned_a),
        "valid": True,
    }


def get_alignment(seq1, seq2, mode='global'):
    """Perform pairwise sequence alignment between two proteins and compute the score using Biopython with the BLOSUM62 matrix."""
    aligner = PairwiseAligner()
    aligner.mode = mode
    try:
        from Bio.Align import substitution_matrices
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5
    except Exception:
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -1.0

    try:
        best_aln = aligner.align(seq1, seq2)[0]
        # Extract aligned sequences with gaps (-) directly from the alignment object
        aligned_h = best_aln[0, :]
        aligned_m = best_aln[1, :]
        return str(best_aln), best_aln.score, aligned_h, aligned_m
    except Exception as error:
        st.warning(f"Sequence alignment failed: {error}")
        return "", 0, "", ""
