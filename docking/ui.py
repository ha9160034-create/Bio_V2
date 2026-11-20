# Docking Streamlit UI Module
import os
import datetime
import hashlib
import json
from pathlib import Path
import streamlit as st
import pandas as pd
import py3Dmol

from docking.models import GridBox
from docking.storage import (
    create_run_directory,
    write_initial_manifest,
    finalize_manifest,
    get_vina_path,
    build_and_save_manifest,
    DISCLAIMER_TEXT,
)
from docking.receptor_prep import prepare_receptor
from docking.ligand_prep import prepare_ligand_from_smiles
from docking.binding_site import (
    define_grid_from_ligand,
    define_grid_from_residues,
    define_grid_whole_protein,
    extract_reference_ligand,
)
from docking.engine import run_vina_multi_seeds, get_vina_version
from docking.validation import (
    calculate_mapped_redocking_rmsd,
    cluster_poses_across_seeds,
    derive_single_docking_status,
    validate_redocking_top_n,
)
from docking.rcsb import fetch_deposited_pdb_text, fetch_biological_assembly_text, fetch_entry_quality
from docking.comparison import (
    run_matched_docking_comparison,
    run_exploratory_cross_docking,
    validate_pocket_chain,
    MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT,
    COMPARISON_LABEL,
    PROTONATION_ASSUMPTION_NOTE,
)
from docking.pair_validation import (
    validate_matched_pair,
    map_healthy_pocket_residues,
    extract_protein_chains,
    MIN_MATCHED_PAIR_IDENTITY,
)
from docking.interactions import (
    calculate_drug_properties,
    analyze_protein_ligand_interactions,
    compare_interactions,
    interpret_score_difference,
    generate_2d_depiction_svg,
    CONTACTS_DISCLAIMER,
)
from docking.openbabel_io import require_openbabel_formats
from docking.pubchem import fetch_drug_from_pubchem
import plotly.graph_objects as go
import numpy as np

# Protocol-specific heuristic bands for the multi-drug panel ΔScore labels.
# These are not binding-affinity thresholds and carry no clinical meaning.
SCORE_DIFF_NEUTRAL_BAND_KCAL = 0.6
SCORE_DIFF_MODERATE_BAND_KCAL = 1.5

def _get_pubchem_cid(smiles: str, state_key: str):
    """Return a PubChem CID only when it belongs to the current SMILES string."""
    info = st.session_state.get(state_key)
    if not isinstance(info, dict):
        return None
    stored_smiles = str(info.get("smiles", "")).strip()
    if stored_smiles and stored_smiles == (smiles or "").strip():
        return info.get("cid")
    return None


def extract_chains_from_pdb(pdb_text: str) -> list[str]:
    chains = set()
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
            ch = line[21].strip()
            if ch:
                chains.add(ch)
    return sorted(list(chains))

def render_docked_complex_3d(receptor_pdb_path: Path, pose_pdbqt_block: str, ref_ligand_pdb_text: str = "",
                              show_surface: bool = False, surface_opacity: float = 0.3, surface_type: str = "MS",
                              pocket_residues: list = None, hbonds: list = None,
                              show_ligand: bool = True, ligand_style: str = "Sticks",
                              grid_box: GridBox = None, show_grid_box: bool = False) -> str:
    view = py3Dmol.view(width="100%", height=480)
    if receptor_pdb_path.exists():
        with open(receptor_pdb_path, "r", encoding="utf-8") as f:
            view.addModel(f.read(), "pdb")
        view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
        if show_surface:
            surf_kind = py3Dmol.MS if surface_type == "MS" else (py3Dmol.VDW if surface_type == "VDW" else py3Dmol.SAS)
            view.addSurface(surf_kind, {"opacity": surface_opacity, "color": "#ECEFF1"}, {"model": 0})
        
        # Highlight binding-pocket residues as sticks
        if pocket_residues:
            for r in pocket_residues:
                r_clean = str(r).strip()
                if r_clean:
                    view.addStyle({"model": 0, "resi": r_clean}, {"stick": {"colorscheme": "orangeCarbon", "radius": 0.25}})

        # Highlight residues participating in polar contacts
        if hbonds:
            for hb in hbonds:
                r_num = str(hb.get("res_num", "")).strip()
                if r_num:
                    view.addStyle({"model": 0, "resi": r_num}, {"stick": {"colorscheme": "magentaCarbon", "radius": 0.32}})

    if ref_ligand_pdb_text and show_ligand:
        view.addModel(ref_ligand_pdb_text, "pdb")
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "grayCarbon", "radius": 0.25}})
    if pose_pdbqt_block and show_ligand:
        view.addModel(pose_pdbqt_block, "pdb")
        if "Spheres" in ligand_style and "Stick" not in ligand_style:
            l_style = {"sphere": {"colorscheme": "greenCarbon", "radius": 0.8}}
        elif "Ball" in ligand_style:
            l_style = {"stick": {"colorscheme": "greenCarbon", "radius": 0.35}, "sphere": {"colorscheme": "greenCarbon", "radius": 0.5}}
        else:
            l_style = {"stick": {"colorscheme": "greenCarbon", "radius": 0.35}}
        view.setStyle({"model": -1}, l_style)
        view.zoomTo({"model": -1})
    else:
        view.zoomTo()

    # Render 3D polygonal wireframe search grid box
    if grid_box and show_grid_box:
        try:
            cx = float(grid_box.center_x)
            cy = float(grid_box.center_y)
            cz = float(grid_box.center_z)
            wx = float(grid_box.size_x)
            hy = float(grid_box.size_y)
            dz = float(grid_box.size_z)

            # 1. Translucent faces to visualize box depth without obscuring the pocket
            view.addBox({
                "center": {"x": cx, "y": cy, "z": cz},
                "dimensions": {"w": wx, "h": hy, "d": dz},
                "color": "#FFD700",
                "opacity": 0.12
            })

            # 2. 3D wireframe cage
            view.addBox({
                "center": {"x": cx, "y": cy, "z": cz},
                "dimensions": {"w": wx, "h": hy, "d": dz},
                "color": "#FFB800",
                "wireframe": True,
                "opacity": 0.85
            })

            # 3. Highlight the 12 boundary edges of the box
            x0, x1 = cx - wx / 2.0, cx + wx / 2.0
            y0, y1 = cy - hy / 2.0, cy + hy / 2.0
            z0, z1 = cz - dz / 2.0, cz + dz / 2.0

            corners = [
                {"x": x0, "y": y0, "z": z0},  # 0
                {"x": x1, "y": y0, "z": z0},  # 1
                {"x": x1, "y": y1, "z": z0},  # 2
                {"x": x0, "y": y1, "z": z0},  # 3
                {"x": x0, "y": y0, "z": z1},  # 4
                {"x": x1, "y": y0, "z": z1},  # 5
                {"x": x1, "y": y1, "z": z1},  # 6
                {"x": x0, "y": y1, "z": z1},  # 7
            ]
            edges = [
                (0, 1), (1, 2), (2, 3), (3, 0),  # Bottom face
                (4, 5), (5, 6), (6, 7), (7, 4),  # Top face
                (0, 4), (1, 5), (2, 6), (3, 7),  # Vertical edges
            ]
            for start_idx, end_idx in edges:
                view.addLine({
                    "start": corners[start_idx],
                    "end": corners[end_idx],
                    "color": "#FF8C00"
                })
        except Exception:
            pass

    return view._make_html()



def sync_target_binding_site_state(current_target_identity: str, session_state_dict: dict) -> bool:
    """
    Syncs and resets target-specific binding site inputs when the target protein changes.
    Returns True if a reset occurred, False otherwise.
    """
    previous_target_identity = session_state_dict.get("dock_previous_target_identity")
    if previous_target_identity != current_target_identity:
        session_state_dict["dock_ref_lig_input"] = (
            "AQ4" if current_target_identity.endswith("1M17") else ""
        )
        session_state_dict["dock_pocket_res_input"] = (
            "790, 858" if current_target_identity.endswith("1M17") else ""
        )
        session_state_dict["dock_pocket_padding"] = 8.0
        session_state_dict["dock_previous_target_identity"] = current_target_identity
        return True
    return False

def find_recovery_pose_index(poses: list, validation_details: dict | None) -> int:
    """
    Finds the index of the pose that actually satisfied the validation gate.

    Full protocol validation points at the top-1 pose; Top-N recovery points at
    the best recovered pose.  Matching is by seed and rank.
    """
    if not poses or not validation_details:
        return 0
    if validation_details.get("protocol_validated"):
        recovered = validation_details.get("top1", {})
    else:
        recovered = validation_details.get("recovery", {})
    recovered_seed = recovered.get("seed")
    recovered_rank = recovered.get("rank")
    if recovered_seed is None or recovered_rank is None:
        return 0
    for index, pose in enumerate(poses):
        if getattr(pose, "seed", None) == recovered_seed and getattr(pose, "rank", None) == recovered_rank:
            return index
    return 0

def compute_input_fingerprint(params: dict) -> str:
    """Compute stable SHA-256 fingerprint for result-defining input parameters."""
    serialized = json.dumps(params, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()

def compute_text_sha256(text: str | None) -> str:
    return hashlib.sha256((text or "").encode("utf-8")).hexdigest()

def lowest_vina_score(poses) -> float | None:
    """Return the lowest Vina score across all calculated poses."""
    scores = [float(pose.score) for pose in poses if getattr(pose, "score", None) is not None]
    return min(scores) if scores else None

def render_docking_tab():
    reason = st.session_state.pop("last_result_invalidation_reason", None)
    if reason:
        st.info(
            f"Previous analysis and docking results were cleared to prevent mixing datasets. "
            f"{reason} Saved run folders remain available as separate provenance records."
        )

    st.header("🔬 Molecular Docking Simulation — Vina")
    vina_ver = get_vina_version()
    st.caption(f"⚙️ Engine: {vina_ver} | Ligand Preparation: RDKit & OpenBabel")
    st.caption(
        "⚠️ Vina scores carry molecular-weight bias; compare only similar ligands and never as experimental ΔG."
    )

    vina_exec = get_vina_path()
    if not vina_exec.exists():
        st.error(f"❌ Vina executable not found at: {vina_exec}")
        return

    try:
        require_openbabel_formats()
    except Exception as exc:
        st.error(f"❌ Open Babel is not ready for docking: {exc}")
        return



    with st.container(border=True):
        tb_col1, tb_col2 = st.columns(2)
        with tb_col1:
            workflow_mode = st.radio(
                "🎯 Analysis Workflow:",
                ["single_docking", "matched_healthy_vs_mutant"],
                format_func=lambda m: "Single-Target Docking" if m == "single_docking" else "Matched WT (Wild-Type Reference) vs. MT (Mutant Variant) Comparison",
                horizontal=True,
                key="docking_workflow_mode"
            )
        with tb_col2:
            ui_mode = st.radio(
                "⚡ Interface Mode:",
                ["⚡ Quick Exploratory Mode", "🔬 Advanced Research Mode"],
                horizontal=True,
                key="docking_interface_mode"
            )
            is_advanced_mode = "Advanced" in ui_mode

        if not is_advanced_mode:
            st.info(
                "⚡ **Quick Exploratory Mode Active:** Streamlined, non-clinical view of the same docking protocol "
                "using standard settings (exhaustiveness = 8, 3 independent seeds)."
            )

    if workflow_mode == "single_docking":
        _render_single_docking_ui(vina_ver, is_advanced_mode=is_advanced_mode)
    else:
        _render_matched_comparison_ui(vina_ver, is_advanced_mode=is_advanced_mode)

def _render_single_docking_ui(vina_ver: str, is_advanced_mode: bool = False):
    svg_2d = ""
    d_props = None

    col_t1, col_t2 = st.columns(2)
    with col_t1:
        with st.container(border=True):
            st.markdown("#### 🏛️ 1️⃣ Receptor Protein & Binding Site Selection")
            target_mode = st.radio(
                "Protein Source:",
                ["custom_pdb", "upload_pdb", "from_tab1"],
                format_func=lambda m: {
                    "custom_pdb": "PDB (Protein Data Bank) ID",
                    "upload_pdb": "Upload PDB (Protein Data Bank) File",
                    "from_tab1": "From Structural Analysis (Tab 1)"
                }[m],
                horizontal=True,
                key="dock_target_source_mode"
            )

            custom_target_id = ""
            uploaded_pdb_file = None
            selected_target_label = ""
            target_key = ""

            if target_mode == "custom_pdb":
                custom_target_id = st.text_input(
                    "PDB (Protein Data Bank) Code:",
                    value=st.session_state.get("dock_custom_pdb_val", ""),
                    placeholder="e.g., 1M17",
                    key="dock_custom_pdb_input"
                ).strip().upper()
                st.session_state["dock_custom_pdb_val"] = custom_target_id
                selected_target_label = custom_target_id if custom_target_id else "CUSTOM"
                target_key = "custom_pdb"
            elif target_mode == "upload_pdb":
                uploaded_pdb_file = st.file_uploader("PDB (Protein Data Bank) File:", type=["pdb"], key="dock_upload_pdb")
                custom_target_id = st.text_input("Protein Name:", value="CUSTOM_RECEPTOR", key="dock_upload_name").strip()
                selected_target_label = custom_target_id
                target_key = "upload_pdb"
            elif target_mode == "from_tab1":
                target_options = []
                if st.session_state.get("h_id"):
                    target_options.append(("WT (Wild-Type Reference): " + st.session_state["h_id"], "h"))
                if st.session_state.get("m_id"):
                    target_options.append(("MT (Mutant Variant): " + st.session_state["m_id"], "m"))
                if target_options:
                    selected_target_label, target_key = st.selectbox(
                        "Select Protein:",
                        target_options,
                        format_func=lambda x: x[0],
                        key="dock_tab1_target_select"
                    )
                else:
                    st.info("No proteins have been loaded in Tab 1 yet.")
                    target_key = "none"
                    selected_target_label = "NONE"

            # Stable identity for the selected target
            if target_mode == "custom_pdb":
                current_target_identity = f"custom:{custom_target_id}"
            elif target_mode == "upload_pdb":
                current_target_identity = (
                    f"upload:{uploaded_pdb_file.name}"
                    if uploaded_pdb_file is not None
                    else "upload:none"
                )
            elif target_mode == "from_tab1":
                current_target_identity = (
                    f"tab1:{target_key}:{st.session_state.get(f'{target_key}_id', '')}"
                )
            else:
                current_target_identity = "custom:none"

            previous_target_identity = st.session_state.get("dock_previous_target_identity")
            was_reset = sync_target_binding_site_state(current_target_identity, st.session_state)
            if was_reset and previous_target_identity is not None:
                st.info("Binding-site settings were reset because the target protein changed.")

            st.markdown("---")
            site_mode = st.radio(
                "Binding-Site Definition:",
                ["co_crystal_ligand", "residue_defined", "blind_docking"],
                format_func=lambda m: (
                    "Co-Crystal Reference Ligand" if m == "co_crystal_ligand"
                    else ("Pocket Residues" if m == "residue_defined"
                    else "🌐 Whole-Protein Blind Docking")
                ),
                horizontal=True,
                key="dock_site_mode_radio"
            )
            keep_crystal_waters = st.checkbox(
                "Keep crystal waters",
                value=False,
                key="dock_keep_waters",
                help=(
                    "Some binding pockets rely on bridging water molecules. When enabled, crystal waters are retained "
                    "in the receptor PDBQT (Protein Data Bank, Partial Charge [Q], & Atom Type [T] Format), "
                    "must be recorded in the run manifest, and must be applied identically to both arms in matched comparisons."
                ),
            )

            ref_ligand_code = ""
            ref_lig_chain = ""
            ref_lig_resseq = ""
            pocket_residues_input = ""
            blind_chain_id = ""
            is_1m17 = (target_mode == "custom_pdb" and custom_target_id == "1M17")

            if site_mode == "co_crystal_ligand":
                col_rf1, col_rf2, col_rf3 = st.columns([1.5, 1.0, 1.0])
                with col_rf1:
                    default_ref_lig = st.session_state.get("dock_ref_lig_input", "") or ("AQ4" if is_1m17 else st.session_state.get("dock_ref_lig_val", ""))
                    ref_ligand_code = st.text_input(
                        "Reference Ligand Code:",
                        value=default_ref_lig,
                        placeholder="e.g., AQ4",
                        key="dock_ref_lig_input"
                    ).strip().upper()
                    if not is_1m17:
                        st.session_state["dock_ref_lig_val"] = ref_ligand_code
                with col_rf2:
                    ref_lig_chain = st.text_input(
                        "Ligand Chain:",
                        value=st.session_state.get("dock_ref_lig_chain", ""),
                        placeholder="e.g., A (optional)",
                        key="dock_ref_lig_chain"
                    ).strip()
                with col_rf3:
                    ref_lig_resseq = st.text_input(
                        "Ligand Residue #:",
                        value=st.session_state.get("dock_ref_lig_resseq", ""),
                        placeholder="e.g., 1001 (optional)",
                        key="dock_ref_lig_resseq"
                    ).strip()
            elif site_mode == "residue_defined":
                default_pocket_res = st.session_state.get("dock_pocket_res_input", "") or ("A:790, A:858" if is_1m17 else st.session_state.get("dock_pocket_res_val", ""))
                pocket_residues_input = st.text_input(
                    "Pocket Residue Numbers (CHAIN:RESIDUE):",
                    value=default_pocket_res,
                    placeholder="e.g., A:790, A:858",
                    key="dock_pocket_res_input"
                )
                if not is_1m17:
                    st.session_state["dock_pocket_res_val"] = pocket_residues_input
            else:
                st.info(
                    "🌐 **Whole-Protein Blind Docking:** A 3D search box enclosing the entire protein structure is generated automatically. "
                    "Raising search exhaustiveness to 16 or higher is recommended for thorough sampling across all surface cavities."
                )
                st.warning(
                    "⚠️ **Low-Confidence Blind Docking:** Results are exploratory and do not substitute for a validated binding site "
                    "(co-crystal ligand or defined pocket residues). The top-scoring pose may occupy a non-functional surface cavity."
                )
                blind_chain_id = st.text_input(
                    "Specific Chain ID (optional; leave blank for whole molecule):",
                    value=st.session_state.get("dock_blind_chain_id", ""),
                    placeholder="e.g., A (optional)",
                    key="dock_blind_chain_id"
                ).strip().upper()

            padding_kwargs = {
                "min_value": 4.0,
                "max_value": 14.0,
                "step": 0.5,
                "format": "%.1f Å",
                "key": "dock_pocket_padding",
            }
            if "dock_pocket_padding" not in st.session_state:
                padding_kwargs["value"] = 8.0
            slider_label = "Box Padding — Å:" if site_mode == "blind_docking" else "Pocket Padding — Å:"
            pocket_padding = st.slider(slider_label, **padding_kwargs)

    with col_t2:
        with st.container(border=True):
            st.markdown("#### 💊 2️⃣ Ligand Properties & Simulation Parameters")
            # Direct live/offline PubChem lookup by generic or trade name
            col_pc1, col_pc2 = st.columns([2.2, 1.0])
            with col_pc1:
                pubchem_query = st.text_input(
                    "🔍 Search Drug in PubChem:",
                    placeholder="e.g., Gefitinib, Aspirin, Imatinib",
                    key="dock_pubchem_search_input"
                )
            with col_pc2:
                st.write("")
                st.write("")
                if st.button("🌐 Fetch Drug", key="btn_pubchem_fetch", use_container_width=True):
                    if pubchem_query.strip():
                        with st.spinner(f"Searching PubChem for '{pubchem_query}'..."):
                            pc_res = fetch_drug_from_pubchem(pubchem_query)
                            if pc_res.get("success"):
                                st.session_state["dock_lig_name_val"] = pc_res["title"]
                                st.session_state["dock_lig_smiles_val"] = pc_res["smiles"]
                                st.session_state["dock_lig_name_input"] = pc_res["title"]
                                st.session_state["dock_lig_smiles_input"] = pc_res["smiles"]
                                st.session_state["dock_pubchem_cid"] = {
                                    "cid": pc_res.get("cid"),
                                    "smiles": pc_res.get("smiles", ""),
                                }
                                st.success(
                                    f"✅ Fetched {pc_res['title']} from {pc_res.get('source')} "
                                    f"(PubChem CID [Compound Identifier]: {pc_res.get('cid')})!"
                                )
                                st.rerun()
                            else:
                                st.error(pc_res.get("error", "PubChem lookup failed."))
                    else:
                        st.warning("Please enter a drug name first.")

            col_dn1, col_dn2 = st.columns([1, 1.6])
            with col_dn1:
                default_lig_name = "Erlotinib" if is_1m17 else st.session_state.get("dock_lig_name_val", "")
                ligand_name = st.text_input(
                    "Drug Name:",
                    value=default_lig_name,
                    placeholder="e.g., Erlotinib",
                    key="dock_lig_name_input"
                )
                if not is_1m17 and ligand_name:
                    st.session_state["dock_lig_name_val"] = ligand_name

            with col_dn2:
                default_smiles = (
                    "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
                    if is_1m17
                    else st.session_state.get("dock_lig_smiles_val", "")
                )
                ligand_smiles = st.text_input(
                    "SMILES (Simplified Molecular-Input Line-Entry System):",
                    value=default_smiles,
                    placeholder="Enter the chemical SMILES string for the drug",
                    key="dock_lig_smiles_input"
                )
                if not is_1m17 and ligand_smiles:
                    st.session_state["dock_lig_smiles_val"] = ligand_smiles

            # 2D chemical depiction and pharmacological properties
            if ligand_smiles.strip():
                tab_struct, tab_lip = st.tabs(["🧪 2D Chemical Structure", "📊 Lipinski Rules — ADME (Absorption, Distribution, Metabolism, and Excretion)"])
                with tab_struct:
                    svg_2d = generate_2d_depiction_svg(ligand_smiles, width=320, height=160)
                    if svg_2d:
                        try:
                            st.image(svg_2d, caption="2D (Two-Dimensional) Chemical Structure Depiction", use_container_width=True)
                        except Exception:
                            st.caption("Could not display 2D chemical depiction.")
                    else:
                        st.caption("Could not render 2D chemical structure.")

                with tab_lip:
                    d_props = calculate_drug_properties(ligand_smiles)
                    if d_props and "error" not in d_props:
                        dp1, dp2, dp3 = st.columns(3)
                        dp1.metric("MW (Molecular Weight)", f"{d_props['molecular_weight']} Da (Daltons / atomic mass units)")
                        dp2.metric("logP (Octanol–Water Partition Coefficient)", f"{d_props['logp']}")
                        dp3.metric("TPSA (Topological Polar Surface Area)", f"{d_props['tpsa']} Å²")
                        dp4, dp5, dp6 = st.columns(3)
                        dp4.metric("HBD (Hydrogen-Bond Donors)", d_props['hbd'])
                        dp5.metric("HBA (Hydrogen-Bond Acceptors)", d_props['hba'])
                        dp6.metric("Rotatable Bonds", d_props['rotatable_bonds'])
                        if d_props['lipinski_pass']:
                            st.caption("✅ Compound complies with Lipinski's Rule of Five for oral drug-likeness.")
                        else:
                            st.caption(f"⚠️ {d_props['drug_likeness']}")
                        if d_props.get("pains_alert"):
                            st.warning(
                                "⚠️ PAINS (Pan-Assay Interference Compounds) Alert: This chemotype may cause assay interference "
                                f"({d_props.get('pains_name', 'PAINS')}). Verify experimentally."
                            )
                        if d_props.get("aggregation_risk"):
                            st.warning(f"⚠️ Heuristic Aggregation Risk: {d_props.get('aggregation_note', '')}")

            st.markdown("---")
            # Search exhaustiveness and random seeds
            if is_advanced_mode:
                st.markdown("##### ⚙️ Advanced Vina Settings:")
                col_cfg1, col_cfg2 = st.columns(2)
                with col_cfg1:
                    exhaustiveness = st.select_slider(
                        "Search Exhaustiveness:",
                        options=[4, 8, 16, 32, 64],
                        value=8,
                        help="Higher values increase binding-pocket sampling thoroughness at the cost of runtime.",
                        key="dock_single_exh"
                    )
                    num_modes = st.slider(
                        "Max Output Modes (num_modes):",
                        min_value=1,
                        max_value=20,
                        value=9,
                        help="Maximum number of binding poses generated per seed.",
                        key="dock_single_num_modes"
                    )
                with col_cfg2:
                    seed_choice = st.selectbox(
                        "Seed Strategy:",
                        ["three_seeds", "single_seed"],
                        format_func=lambda x: "3 Seeds [42, 101, 2024] (Statistical Reproducibility)" if x == "three_seeds" else "Single Seed [42] (Quick Test)",
                        key="dock_single_seeds"
                    )
                    energy_range = st.slider(
                        "Pose Energy Range — kcal/mol:",
                        min_value=1.0,
                        max_value=5.0,
                        value=3.0,
                        step=0.5,
                        help="Maximum energy difference above the best pose to retain a binding mode.",
                        key="dock_single_energy_range"
                    )
            else:
                exhaustiveness = 8
                seed_choice = "three_seeds"
                num_modes = 9
                energy_range = 3.0
                st.caption("⚡ **Standard Quick Settings:** Exhaustiveness = 8 | 3 independent random seeds [42, 101, 2024] for reproducibility | 9 modes.")

    single_fp_payload = {
        "workflow": "single_docking",
        "target_identity": current_target_identity,
        "site_mode": site_mode,
        "ref_ligand_code": ref_ligand_code if site_mode == "co_crystal_ligand" else "",
        "ref_lig_chain": ref_lig_chain if site_mode == "co_crystal_ligand" else "",
        "ref_lig_resseq": ref_lig_resseq if site_mode == "co_crystal_ligand" else "",
        "pocket_residues": pocket_residues_input if site_mode == "residue_defined" else "",
        "blind_chain_id": blind_chain_id if site_mode == "blind_docking" else "",
        "pocket_padding": pocket_padding,
        "ligand_name": ligand_name,
        "ligand_smiles": ligand_smiles,
        "exhaustiveness": exhaustiveness,
        "seed_choice": seed_choice,
        "num_modes": num_modes,
        "energy_range": energy_range,
        "vina_version": vina_ver,
    }
    current_single_fp = compute_input_fingerprint(single_fp_payload)
    prev_single_fp = st.session_state.get("dock_single_input_fingerprint")
    if prev_single_fp and prev_single_fp != current_single_fp:
        st.session_state.pop("current_docking_result", None)
    st.session_state["dock_single_input_fingerprint"] = current_single_fp

    if st.button("🚀 Run Docking", type="primary"):
        st.session_state.pop("current_docking_result", None)
        if site_mode == "co_crystal_ligand" and not ref_ligand_code:
            st.error("Please enter a reference ligand code to define the binding site.")
            return
        if site_mode == "residue_defined" and not pocket_residues_input.strip():
            st.error("Please enter binding-site residue numbers.")
            return
        if not ligand_smiles.strip():
            st.error("Please enter a SMILES (Simplified Molecular-Input Line-Entry System) string for the ligand.")
            return

        seeds_to_run = [42, 101, 2024] if seed_choice == "three_seeds" else [42]
        created_at_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()

        # 1. Coordinate retrieval and exact provenance
        if target_key == "demo_1M17":
            target_id = "1M17"
            assembly_text = fetch_biological_assembly_text("1M17", 1)
            if assembly_text:
                pdb_data = assembly_text
                coord_source = "rcsb_biological_assembly"
                assembly_id = 1
                source_url = "https://files.rcsb.org/download/1M17.pdb1.gz"
            else:
                pdb_data = fetch_deposited_pdb_text("1M17")
                coord_source = "deposited"
                assembly_id = None
                source_url = "https://files.rcsb.org/download/1M17.pdb"
        elif target_key == "custom_pdb":
            if not custom_target_id:
                st.error("Please enter a valid PDB (Protein Data Bank) code.")
                return
            target_id = custom_target_id
            with st.spinner(f"Fetching protein structure {target_id} from RCSB PDB (Research Collaboratory for Structural Bioinformatics Protein Data Bank)..."):
                assembly_text = fetch_biological_assembly_text(target_id, 1)
                if assembly_text:
                    pdb_data = assembly_text
                    coord_source = "rcsb_biological_assembly"
                    assembly_id = 1
                    source_url = f"https://files.rcsb.org/download/{target_id}.pdb1.gz"
                else:
                    pdb_data = fetch_deposited_pdb_text(target_id)
                    coord_source = "deposited"
                    assembly_id = None
                    source_url = f"https://files.rcsb.org/download/{target_id}.pdb"
            if not pdb_data or len(pdb_data.strip()) == 0:
                st.error(f"❌ Could not fetch PDB (Protein Data Bank) entry {target_id}. Please check the code or your internet connection.")
                return
        elif target_key == "upload_pdb":
            if not uploaded_pdb_file:
                st.error("Please upload a PDB (Protein Data Bank) file first.")
                return
            target_id = custom_target_id if custom_target_id else "UPLOADED_PROTEIN"
            pdb_data = uploaded_pdb_file.getvalue().decode("utf-8", errors="ignore")
            coord_source = "user_upload"
            assembly_id = None
            source_url = uploaded_pdb_file.name
        else:
            target_id = st.session_state.get(f"{target_key}_id", "PROTEIN")
            if st.session_state.get(f"{target_key}_assembly_pdb"):
                pdb_data = st.session_state[f"{target_key}_assembly_pdb"]
                coord_source = "rcsb_biological_assembly"
                assembly_id = 1
                source_url = f"https://files.rcsb.org/download/{target_id}.pdb1.gz"
            elif st.session_state.get(f"{target_key}_pdb"):
                pdb_data = st.session_state[f"{target_key}_pdb"]
                coord_source = "deposited"
                assembly_id = None
                source_url = f"https://files.rcsb.org/download/{target_id}.pdb"
            else:
                st.error("Coordinates for the selected protein were not found.")
                return

        chains = extract_chains_from_pdb(pdb_data)
        quality_meta = (
            fetch_entry_quality(
                target_id,
                ligand_code=ref_ligand_code if site_mode == "co_crystal_ligand" else "",
                include_residue_rscc=True,
            )
            if coord_source in {"deposited", "rcsb_biological_assembly"}
            else {}
        )
        structure_advisories = []
        low_residue_rscc = []
        if quality_meta:
            quality_bits = [
                f"Resolution: {quality_meta.get('resolution_angstrom') or 'unknown'} Å",
                f"R_free: {quality_meta.get('r_free') if quality_meta.get('r_free') is not None else 'unknown'}",
                f"Method: {quality_meta.get('experimental_method') or 'unknown'}",
            ]
            if quality_meta.get("ligand_rscc") is not None:
                quality_bits.append(
                    f"Ref-ligand RSCC: {quality_meta['ligand_rscc']:.3f} "
                    f"({quality_meta.get('ligand_code', '')})"
                )
            residue_rscc = quality_meta.get("residue_rscc_summary") or {}
            if residue_rscc.get("n_residues"):
                quality_bits.append(
                    f"Protein RSCC median: {residue_rscc['median_rscc']:.3f} "
                    f"(<0.8: {residue_rscc['n_below_0_8']}/{residue_rscc['n_residues']})"
                )
            st.caption(" | ".join(quality_bits))
            if (
                quality_meta.get("resolution_angstrom") is not None
                and quality_meta["resolution_angstrom"] > 2.5
            ):
                structure_advisories.append(
                    "⚠️ Low structural resolution (> 2.5 Å); binding-pocket geometry may be uncertain. "
                    "This is a structure-quality advisory, not an automatic rejection."
                )
            if quality_meta.get("r_free") is not None and quality_meta["r_free"] >= 0.25:
                structure_advisories.append(
                    f"⚠️ R_free = {quality_meta['r_free']:.3f} (≥ 0.25); model quality is below the "
                    "benchmark guideline (R_free < 0.25). Binding-pocket geometry and ligand "
                    "placement may be less reliable. This is an advisory, not an automatic rejection."
                )
            if (
                quality_meta.get("ligand_rscc") is not None
                and quality_meta["ligand_rscc"] < 0.8
            ):
                structure_advisories.append(
                    "⚠️ Reference-ligand RSCC < 0.8; electron density is weak "
                    "and the crystal pose may be uncertain."
                )
            if (
                residue_rscc.get("median_rscc") is not None
                and (
                    residue_rscc["median_rscc"] < 0.8
                    or (residue_rscc.get("fraction_below_0_8") or 0) > 0.1
                )
            ):
                structure_advisories.append(
                    "⚠️ Protein-residue RSCC is low (median < 0.8 or >10% of residues < 0.8); "
                    "review local density quality before relying on the pocket."
                )
            low_residue_rscc = quality_meta.get("residue_rscc_low") or []
        run_id, run_dir = create_run_directory(f"dock_{target_id}")

        # Save immutable source PDB
        source_pdb_path = run_dir / "receptor_source.pdb"
        with open(source_pdb_path, "w", encoding="utf-8") as f:
            f.write(pdb_data)

        target_meta = {
            "pdb_id": target_id,
            "coordinate_source": coord_source,
            "assembly_id": assembly_id,
            "source_url": source_url,
            "chains": chains,
            "structure_quality": quality_meta,
            "created_at_utc": created_at_utc,
            "mutation_label": st.session_state.get(f"{target_key}_mutation_label", "")
        }
        binding_site_meta = {
            "method": site_mode,
            "reference_ligand": ref_ligand_code if site_mode == "co_crystal_ligand" else None,
            "reference_ligand_chain": ref_lig_chain if site_mode == "co_crystal_ligand" and ref_lig_chain else None,
            "reference_ligand_resseq": ref_lig_resseq if site_mode == "co_crystal_ligand" and ref_lig_resseq else None,
            "residues": pocket_residues_input if site_mode == "residue_defined" else None,
            "chain_id": blind_chain_id if site_mode == "blind_docking" and blind_chain_id else None,
            "padding": pocket_padding
        }
        ligand_meta = {
            "name": ligand_name,
            "input_smiles": ligand_smiles,
            "protonation_target_ph": 7.4,
        }
        engine_meta = {
            "name": "AutoDock Vina",
            "version": vina_ver,
            "seeds": seeds_to_run,
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "energy_range": energy_range
        }

        # 2. Binding site definition (Prohibits origin fallback)
        ref_ligand_pdb_text = ""
        try:
            if site_mode == "co_crystal_ligand":
                parsed_resseq = int(ref_lig_resseq) if ref_lig_resseq and ref_lig_resseq.isdigit() else None
                grid, ref_coords, ref_ligand_pdb_text = define_grid_from_ligand(
                    pdb_data,
                    ref_ligand_code,
                    padding=pocket_padding,
                    chain_id=ref_lig_chain or None,
                    resseq=parsed_resseq
                )
                with open(run_dir / "reference_ligand.pdb", "w", encoding="utf-8") as f:
                    f.write(ref_ligand_pdb_text + "\n")
            elif site_mode == "residue_defined":
                res_specs = [r.strip() for r in pocket_residues_input.split(",") if r.strip()]
                grid = define_grid_from_residues(pdb_data, res_specs, padding=pocket_padding)
            elif site_mode == "blind_docking":
                grid = define_grid_whole_protein(
                    pdb_data,
                    padding=pocket_padding,
                    chain_id=blind_chain_id or None
                )
        except Exception as e:
            manifest_fail = {
                "schema_version": 1,
                "run_id": run_id,
                "status": "run failed",
                "created_at_utc": created_at_utc,
                "error": str(e),
                "target": target_meta,
                "binding_site": binding_site_meta,
                "ligand": ligand_meta,
                "engine": engine_meta,
                "warnings": [f"Binding site error: {e}"]
            }
            finalize_manifest(run_dir, manifest_fail)
            st.error(f"Docking was not run: define a validated binding site. Error: {e}")
            return

        grid_meta = {
            "center": [grid.center_x, grid.center_y, grid.center_z],
            "size": [grid.size_x, grid.size_y, grid.size_z]
        }

        grid_box_warnings = []
        oversized_axes = [
            f"{axis} = {value:.1f} Å"
            for axis, value in zip("xyz", (grid.size_x, grid.size_y, grid.size_z))
            if value > 30.0
        ]
        if oversized_axes:
            grid_warning = (
                "⚠️ Search box exceeds the 30 Å guideline on "
                + ", ".join(oversized_axes)
                + ". Prefer the smallest box that contains the pocket; reduce padding or "
                "increase exhaustiveness."
            )
            structure_advisories.append(grid_warning)
            grid_box_warnings.append(grid_warning)

        if structure_advisories:
            with st.expander("⚠️ Structure Quality Advisories", expanded=False):
                for adv in structure_advisories:
                    st.warning(adv)
                if low_residue_rscc:
                    with st.expander("🔎 Protein Residues with RSCC < 0.8"):
                        st.dataframe(
                            low_residue_rscc,
                            use_container_width=True,
                            hide_index=True,
                        )

        # 3. Write initial manifest before running executables
        initial_data = {
            "schema_version": 1,
            "run_id": run_id,
            "status": "running",
            "created_at_utc": created_at_utc,
            "target": target_meta,
            "binding_site": binding_site_meta,
            "ligand": ligand_meta,
            "grid": grid_meta,
            "engine": engine_meta
        }
        write_initial_manifest(run_dir, initial_data)

        # 4. Receptor & Ligand Preparation
        with st.spinner("Preparing receptor and ligand and running Vina..."):
            try:
                clean_pdb, rec_pdbqt, prep_rep = prepare_receptor(
                    pdb_data,
                    run_dir,
                    keep_waters=keep_crystal_waters,
                    reference_ligand_resname=ref_ligand_code if site_mode == "co_crystal_ligand" else None
                )
                sdf_path, lig_pdbqt, lig_rep = prepare_ligand_from_smiles(
                    ligand_smiles,
                    ligand_name,
                    run_dir,
                    ph=7.4,
                    pubchem_cid=_get_pubchem_cid(ligand_smiles, "dock_pubchem_cid"),
                )
            except Exception as e:
                build_and_save_manifest(
                    run_dir=run_dir,
                    run_id=run_id,
                    status="run failed",
                    target_meta=target_meta,
                    binding_site_meta=binding_site_meta,
                    ligand_meta=ligand_meta,
                    prep_meta={},
                    engine_meta=engine_meta,
                    grid_meta=grid_meta,
                    results_meta={},
                    warnings=[f"Preparation failure: {e}"]
                )
                st.error(f"Input preparation failed: {e}")
                return

            # 5. Execute Multi-Seed Vina runs
            try:
                poses, runs_meta = run_vina_multi_seeds(
                    receptor_pdbqt=rec_pdbqt,
                    ligand_pdbqt=lig_pdbqt,
                    grid=grid,
                    output_dir=run_dir,
                    seeds=seeds_to_run,
                    exhaustiveness=exhaustiveness,
                    num_modes=num_modes,
                    energy_range=energy_range
                )
            except Exception as e:
                build_and_save_manifest(
                    run_dir=run_dir,
                    run_id=run_id,
                    status="run failed",
                    target_meta=target_meta,
                    binding_site_meta=binding_site_meta,
                    ligand_meta=lig_rep,
                    prep_meta=prep_rep,
                    engine_meta=engine_meta,
                    grid_meta=grid_meta,
                    results_meta={},
                    warnings=[f"Engine failure: {e}"]
                )
                st.error(f"Vina execution failed: {e}")
                return

            # 6. Cluster poses across seeds
            clusters = cluster_poses_across_seeds(poses, rmsd_threshold=2.0)

            # 7. Redocking validation (Top-N validation per REDOCKING_VALIDATION_TOP_N_REPAIR)
            VALIDATION_TOP_N = 10
            MIN_VALIDATION_SEEDS = 3
            MIN_CLUSTER_SEED_SUPPORT = 2
            RMSD_CUTOFF_ANGSTROM = 2.0

            redocking_rmsd = None
            validation_status = "exploratory—no reference validation"
            warnings = list(grid_box_warnings)
            validated_pose_info = {}
            mapping_details = {}
            validation_details = None

            if len(seeds_to_run) < MIN_VALIDATION_SEEDS:
                validation_status = "exploratory—insufficient independent seeds"
                warnings.append("Single-seed/insufficient seeds: pose reproducibility across independent runs cannot be verified.")

            if site_mode == "co_crystal_ligand" and ref_ligand_pdb_text:
                try:
                    validation_details = validate_redocking_top_n(
                        ref_ligand_pdb_text,
                        poses,
                        clusters,
                        top_n=VALIDATION_TOP_N,
                        rmsd_cutoff=RMSD_CUTOFF_ANGSTROM,
                        min_seed_support=MIN_CLUSTER_SEED_SUPPORT,
                    )

                    top1_info = validation_details["top1"]
                    rec_info = validation_details["recovery"]

                    redocking_rmsd = top1_info["rmsd_angstrom"]
                    mapping_details = top1_info["mapping_details"]

                    validation_status, gate_warnings = derive_single_docking_status(
                        validation_details,
                        len(seeds_to_run),
                        min_validation_seeds=MIN_VALIDATION_SEEDS,
                        top_n=VALIDATION_TOP_N,
                        rmsd_cutoff=RMSD_CUTOFF_ANGSTROM,
                    )
                    warnings.extend(gate_warnings)

                    chosen_info = None
                    if validation_details.get("protocol_validated"):
                        chosen_info = top1_info
                    elif (
                        validation_details.get("recovery_rmsd_pass")
                        and validation_details.get("recovery_cluster_pass")
                    ):
                        chosen_info = rec_info
                    if chosen_info is not None:
                        validated_pose_info = {
                            "seed": chosen_info["seed"],
                            "rank": chosen_info["rank"],
                            "cluster_id": chosen_info["cluster_id"],
                            "output_file": f"out_seed_{chosen_info['seed']}.pdbqt",
                            "score": chosen_info["score"],
                            "rmsd_angstrom": chosen_info["rmsd_angstrom"],
                            "cluster_seed_count": chosen_info["cluster_seed_count"],
                            "supporting_seeds": chosen_info["supporting_seeds"]
                        }
                except Exception as ex:
                    validation_status = "validation unavailable"
                    warnings.append(f"Redocking calculation unavailable: {ex}")

            # Calculate global top_score across every pose from every seed
            global_top_score = min((p.score for p in poses), default=None)

            redocking_manifest_data = {
                "reference_ligand_file": "reference_ligand.pdb" if ref_ligand_pdb_text else None,
                "validated_pose": validated_pose_info,
                "rmsd_angstrom": redocking_rmsd,
                "mapping_details": mapping_details
            }
            if validation_details is not None:
                redocking_manifest_data.update({
                    "criterion": {
                        "top_n": VALIDATION_TOP_N,
                        "rmsd_cutoff_angstrom": RMSD_CUTOFF_ANGSTROM,
                        "minimum_independent_seeds": MIN_VALIDATION_SEEDS,
                        "minimum_cluster_seed_support": MIN_CLUSTER_SEED_SUPPORT
                    },
                    "top1_scored_pose": {
                        "seed": validation_details["top1"]["seed"],
                        "rank": validation_details["top1"]["rank"],
                        "score": validation_details["top1"]["score"],
                        "cluster_id": validation_details["top1"]["cluster_id"],
                        "output_file": f"out_seed_{validation_details['top1']['seed']}.pdbqt",
                        "rmsd_angstrom": validation_details["top1"]["rmsd_angstrom"],
                        "mapping_details": validation_details["top1"]["mapping_details"]
                    },
                    "best_recovered_pose_in_top_n": {
                        "seed": validation_details["recovery"]["seed"],
                        "rank": validation_details["recovery"]["rank"],
                        "score": validation_details["recovery"]["score"],
                        "cluster_id": validation_details["recovery"]["cluster_id"],
                        "output_file": f"out_seed_{validation_details['recovery']['seed']}.pdbqt",
                        "rmsd_angstrom": validation_details["recovery"]["rmsd_angstrom"],
                        "cluster_seed_count": validation_details["recovery"]["cluster_seed_count"],
                        "supporting_seeds": validation_details["recovery"]["supporting_seeds"],
                        "mapping_details": validation_details["recovery"]["mapping_details"]
                    },
                    "top1_rmsd_pass": validation_details["top1_rmsd_pass"],
                    "top1_cluster_pass": validation_details["top1_cluster_pass"],
                    "recovery_rmsd_pass": validation_details["recovery_rmsd_pass"],
                    "recovery_cluster_pass": validation_details["recovery_cluster_pass"],
                    "protocol_validated": validation_details["protocol_validated"],
                    "validation_status": validation_status,
                })

            results_meta = {
                "top_score": global_top_score,
                "validation_status": validation_status,
                "per_seed": runs_meta,
                "clusters": [
                    {
                        "cluster_id": cl.cluster_id,
                        "supporting_seeds": cl.supporting_seeds,
                        "seed_count": cl.seed_count,
                        "top_score": cl.top_score,
                        "median_score": cl.median_score,
                        "representative_pose_source": cl.representative_pose_source
                    }
                    for cl in clusters
                ],
                "redocking": redocking_manifest_data
            }

            # 8. Finalize Manifest atomically
            final_manifest_path = build_and_save_manifest(
                run_dir=run_dir,
                run_id=run_id,
                status=validation_status,
                target_meta=target_meta,
                binding_site_meta=binding_site_meta,
                ligand_meta=lig_rep,
                prep_meta=prep_rep,
                engine_meta=engine_meta,
                grid_meta=grid_meta,
                results_meta=results_meta,
                warnings=warnings,
                input_fingerprint=current_single_fp
            )

            st.session_state["current_docking_result"] = {
                "run_id": run_id,
                "target_id": target_id,
                "ligand_name": ligand_name,
                "clean_pdb": clean_pdb,
                "prep_report": prep_rep,
                "poses": poses,
                "clusters": clusters,
                "grid": grid,
                "ref_ligand_pdb": ref_ligand_pdb_text,
                "redocking_rmsd": redocking_rmsd,
                "validation_status": validation_status,
                "validation_details": validation_details,
                "warnings": warnings,
                "run_dir": run_dir,
                "input_fingerprint": current_single_fp
            }
            st.success(f"✅ Docking completed successfully! ({run_id})")

    # Display results
    res = st.session_state.get("current_docking_result")
    if res:
        st.markdown("---")
        st.subheader(f"📊 Docking Results: {res['ligand_name']} ⚯ {res['target_id']}")

        val_stat = res["validation_status"]
        val_det = res.get("validation_details")
        best_score = lowest_vina_score(res.get("poses", []))
        best_score_display = f"{best_score:.2f} kcal/mol" if best_score is not None else "—"
        total_poses_count = len(res.get("poses", []))
        clusters_count = len(res.get("clusters", []))

        # Top KPI Metric Cards
        with st.container(border=True):
            kpi_col1, kpi_col2, kpi_col3, kpi_col4 = st.columns(4)
            if val_stat in {"validation passed", "protocol validated on WT reference"}:
                gate_badge = "🟢 Validation Passed"
            elif val_stat == "validation passed (Top-10 recovery)":
                gate_badge = "🟡 Top-10 Recovery Only (Not Top-1)"
            elif "failed" in val_stat:
                gate_badge = "🔴 Validation Failed"
            elif "insufficient" in val_stat:
                gate_badge = "🟡 Exploratory (Insufficient Seeds)"
            else:
                gate_badge = "⚪ Unavailable (No Reference Ligand)"

            kpi_col1.metric("🎯 Validation Status", gate_badge)
            kpi_col2.metric("Lowest Vina Score", best_score_display)
            kpi_col3.metric("🔢 Calculated Poses", f"{total_poses_count} poses")
            kpi_col4.metric("👥 Pose Clusters", f"{clusters_count} clusters")

        # Exploratory Summary for Quick Mode
        if not is_advanced_mode:
            with st.container(border=True):
                st.markdown("#### Exploratory Docking Summary")
                if best_score is None:
                    st.markdown("- No docked pose is available.")
                else:
                    st.markdown(
                        f"- **Lowest Vina score:** `{best_score:.2f} kcal/mol`"
                    )
                if d_props and "lipinski_pass" in d_props:
                    st.markdown(
                        f"- **Rule-of-Five:** `{d_props['lipinski_violations']}` violations"
                    )

        # Status badge & warnings
        if val_stat == "validation passed (Top-10 recovery)":
            rec_r = val_det['recovery']['rmsd_angstrom'] if val_det else res['redocking_rmsd']
            st.warning(
                f"🟡 Validation Status: {val_stat} — Best Top-10 recovery RMSD = {rec_r} Å <= 2.0 Å. "
                "Top-1 did not reproduce the crystal pose; this is not full protocol validation."
            )
        elif val_stat == "validation passed":
            st.success(f"🎯 Validation Status: validation passed — RMSD = {res['redocking_rmsd']} Å <= 2.0 Å")
        elif val_stat == "validation failed":
            st.error("❌ Validation Status: validation failed")
        elif val_stat == "validation unavailable":
            st.warning("⚠️ Validation Status: validation unavailable (atom count or chemical identity does not match the reference ligand).")
        else:
            st.info(f"ℹ️ Status: {val_stat}")

        # Exact user-facing metrics: shown directly in Advanced Mode, folded in Quick Mode
        if val_det:
            t1 = val_det["top1"]
            rec = val_det["recovery"]
            val_summary_md = (
                f"- **Top-1 Vina pose RMSD:** `{t1['rmsd_angstrom']:.3f} Å` — Score: `{t1['score']} kcal/mol`\n"
                f"- **Best Top-10 recovery RMSD:** `{rec['rmsd_angstrom']:.3f} Å` — Seed `{rec['seed']}`, Rank `{rec['rank']}`, Score: `{rec['score']} kcal/mol`\n"
                f"- **Recovery cluster support:** `{rec['cluster_seed_count']}` independent seeds `{rec['supporting_seeds']}`"
            )
            if is_advanced_mode:
                st.markdown(val_summary_md)
            else:
                with st.expander("🔬 Pose Validation Details"):
                    st.markdown(val_summary_md)

        for w in res.get("warnings", []):
            clean_w = w.lstrip("⚠️ ").strip()
            if "Search box exceeds the 30 Å guideline" in clean_w:
                continue
            st.warning(f"⚠️ {clean_w}")

        # Organized Result Tabs
        tab_3d, tab_poses = st.tabs([
            "🔮 3D (Three-Dimensional) View & Contacts",
            "📋 Poses & Clusters Table"
        ])

        with tab_3d:
            # 3D Toolbar controls
            with st.container(border=True):
                c_row1, c_row2, c_row3, c_row4 = st.columns([1.5, 1.2, 1.0, 1.0])
                with c_row1:
                    default_pose_idx = find_recovery_pose_index(res["poses"], res.get("validation_details"))
                    sel_pose_idx = st.selectbox(
                        "Select Pose to Display:",
                        range(len(res["poses"])),
                        index=default_pose_idx,
                        key=f"dock_pose_selector_{res['run_id']}",
                        format_func=lambda i: (
                            f"Pose #{res['poses'][i].rank} [Seed {res['poses'][i].seed}] "
                            f"({res['poses'][i].score} kcal/mol)"
                        ),
                    )
                with c_row2:
                    dock_show_ligand = st.checkbox("Show Ligand", value=True, key="dock_single_show_ligand")
                    if dock_show_ligand:
                        dock_lig_style = st.selectbox("Ligand Style:", ["Sticks", "Ball & Stick", "Spheres"], key="dock_single_lig_style")
                    else:
                        dock_lig_style = "Sticks"
                with c_row3:
                    dock_show_grid = st.checkbox("📦 Search Grid Box", value=False, key="dock_single_show_grid")
                    dock_show_surface = st.checkbox("Show Surface", value=False, key="dock_single_show_surface")
                with c_row4:
                    if dock_show_surface:
                        dock_surf_type = st.selectbox("Surface Type:", ["MS (Molecular Surface)", "SAS (Solvent-Accessible Surface)", "VDW (van der Waals Surface)"], key="dock_single_surf_type")
                        dock_surf_opacity = st.slider("Surface Opacity:", 0.0, 1.0, 0.3, key="dock_single_surf_opacity")
                    else:
                        dock_surf_type = "MS"
                        dock_surf_opacity = 0.3

            validation_details = res.get("validation_details")
            if validation_details and "recovery" in validation_details:
                if validation_details.get("protocol_validated"):
                    validated_pose = validation_details.get("top1", {})
                    validated_label = "validated top-1 pose"
                else:
                    validated_pose = validation_details["recovery"]
                    validated_label = "validated Top-10 recovery pose"
                sel_p = res["poses"][sel_pose_idx]
                if getattr(sel_p, "seed", None) == validated_pose.get("seed") and getattr(sel_p, "rank", None) == validated_pose.get("rank"):
                    st.caption(
                        f"🎯 Displayed pose: {validated_label} — "
                        f"RMSD = {validated_pose['rmsd_angstrom']:.3f} Å."
                    )
                else:
                    st.caption("ℹ️ User-selected pose displayed; not necessarily the validated pose.")

            # Analyze interactions for selected pose
            active_pose = res["poses"][sel_pose_idx]
            try:
                with open(res["clean_pdb"], "r", encoding="utf-8") as f:
                    rec_clean_text = f.read()
                single_inter = analyze_protein_ligand_interactions(
                    receptor_pdb_text=rec_clean_text,
                    ligand_pdbqt_block=active_pose.pdbqt_block
                )
            except Exception:
                single_inter = {
                    "polar_contact_candidates": [],
                    "carbon_contact_candidates": [],
                    "interacting_residues": [],
                    "polar_contact_count": 0,
                    "carbon_contact_count": 0,
                    "hbonds": [],
                    "hydrophobic": [],
                    "hbond_count": 0,
                    "hydrophobic_count": 0,
                }

            # Layout: 60% 3D Viewer, 40% Interactions Details
            col_v3d, col_vint = st.columns([1.5, 1.0])
            with col_v3d:
                st.caption("🟢 Green: Calculated Ligand Pose | ⚪ Gray: Crystal Reference Ligand | 🟣 Magenta: Polar-Contact Residues")
                active_pose_pdbqt = active_pose.pdbqt_block
                polar_candidates = single_inter.get("polar_contact_candidates") or single_inter.get("hbonds", [])
                view_html = render_docked_complex_3d(
                    res["clean_pdb"],
                    active_pose_pdbqt,
                    res["ref_ligand_pdb"],
                    show_surface=dock_show_surface,
                    surface_opacity=dock_surf_opacity,
                    surface_type=dock_surf_type.split()[0],
                    hbonds=polar_candidates,
                    show_ligand=dock_show_ligand,
                    ligand_style=dock_lig_style,
                    grid_box=res.get("grid"),
                    show_grid_box=dock_show_grid
                )
                st.iframe(view_html, height=520)

            with col_vint:
                with st.container(border=True):
                    st.markdown("##### 🔗 Exploratory Pose Contacts")
                    st.info(CONTACTS_DISCLAIMER)
                    si_col1, si_col2 = st.columns(2)
                    p_count = single_inter.get("polar_contact_count", single_inter.get("hbond_count", 0))
                    c_count = single_inter.get("carbon_contact_count", single_inter.get("hydrophobic_count", 0))
                    si_col1.metric("Polar Contact Candidates", p_count)
                    si_col2.metric("Carbon Contact Candidates", c_count)

                    if polar_candidates:
                        st.markdown("###### 💧 Polar Contact Candidates (<= 3.5 Å):")
                        df_hb = pd.DataFrame(polar_candidates).rename(columns={
                            "residue": "Residue",
                            "chain": "Chain",
                            "rec_atom": "Receptor Atom",
                            "lig_atom": "Ligand Atom",
                            "distance_angstrom": "Distance — Å"
                        })[["Residue", "Chain", "Receptor Atom", "Ligand Atom", "Distance — Å"]]
                        st.dataframe(df_hb, use_container_width=True, hide_index=True, height=160)
                    else:
                        st.info("No direct polar contact candidates within <= 3.5 Å for this pose.")

                    if single_inter.get("interacting_residues"):
                        st.markdown("###### 🧱 Surrounding Pocket Residues (<= 4.0 Å):")
                        df_ir = pd.DataFrame(single_inter["interacting_residues"]).rename(columns={
                            "res_name": "Residue",
                            "res_num": "Number",
                            "chain": "Chain",
                            "min_distance": "Min Distance — Å"
                        })[["Residue", "Number", "Chain", "Min Distance — Å"]]
                        st.dataframe(df_ir, use_container_width=True, hide_index=True, height=160)

        with tab_poses:
            poses_data = [
                {
                    "Rank": p.rank,
                    "Vina Score — kcal/mol": p.score,
                    "Seed": p.seed,
                    "Cluster": p.cluster_id,
                    "RMSD l.b.": p.rmsd_lb,
                    "RMSD u.b.": p.rmsd_ub,
                    "Reference RMSD — Å": f"{p.rmsd_to_reference:.3f}" if p.rmsd_to_reference is not None else "—"
                }
                for p in res["poses"]
            ]
            df_poses = pd.DataFrame(poses_data)

            st.markdown("#### 📋 Calculated Poses Table")
            st.dataframe(df_poses, use_container_width=True, hide_index=True)

            if res.get("clusters"):
                st.markdown("#### 👥 Pose Clusters")
                cl_data = [
                    {
                        "Cluster": cl.cluster_id,
                        "Supporting Seeds Count": cl.seed_count,
                        "Seeds": str(cl.supporting_seeds),
                        "Top Score — kcal/mol": cl.top_score,
                        "Median Score — kcal/mol": cl.median_score
                    }
                    for cl in res["clusters"]
                ]
                st.dataframe(pd.DataFrame(cl_data), use_container_width=True, hide_index=True)



def _render_matched_comparison_ui(vina_ver: str, is_advanced_mode: bool = False):
    st.subheader("⚖️ Comparative Docking Analysis — WT (Wild-Type Reference) vs. MT (Mutant Variant)")

    h_active_pdb = st.session_state.get("h_assembly_pdb") or st.session_state.get("h_pdb") or ""
    m_active_pdb = st.session_state.get("m_assembly_pdb") or st.session_state.get("m_pdb") or ""

    if not h_active_pdb or not m_active_pdb:
        st.warning("Please load both WT (Wild-Type Reference) and MT (Mutant Variant) structures in Tab 1 (Structural Analysis) first.")
        return

    # Extract available protein chains
    h_chains = extract_protein_chains(h_active_pdb)
    m_chains = extract_protein_chains(m_active_pdb)

    h_chain_default = st.session_state.get("h_selected_chain")
    if h_chain_default not in h_chains and h_chains:
        h_chain_default = h_chains[0]
    m_chain_default = st.session_state.get("m_selected_chain")
    if m_chain_default not in m_chains and m_chains:
        m_chain_default = m_chains[0]

    # 1. Workflow Selection
    workflow_mode = st.radio(
        "Comparison Workflow:",
        ["matched_variant", "exploratory_cross"],
        format_func=lambda w: (
            "⚖️ Matched WT (Wild-Type Reference) vs. MT (Mutant Variant) Comparison [Requires >= 90% Sequence Identity]"
            if w == "matched_variant"
            else "🔀 Exploratory Cross-Protein Docking [Independent Proteins/Pockets; No ΔScore (Docking Score Difference)]"
        ),
        horizontal=True,
        key="comp_workflow_choice",
    )

    # Optional documented example preset (never hardcoded as default)
    with st.expander("📋 Load Documented Example (Optional)", expanded=False):
        st.info(
            "This is a documented illustrative example only for the EGFR (Epidermal Growth Factor Receptor) Kinase Domain "
            "with Erlotinib and the T790M mutation. It is not a generic default for other proteins."
        )
        if st.button("Load EGFR (Epidermal Growth Factor Receptor) Documented Example", key="btn_load_egfr_preset"):
            st.session_state["comp_pocket_val"] = "A:766, A:769, A:834"
            st.session_state["comp_domain_val"] = "Kinase Domain"
            st.session_state["comp_lig_name_val"] = "Erlotinib"
            st.session_state["comp_smiles_val"] = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
            st.session_state["comp_pocket_input"] = "A:766, A:769, A:834"
            st.session_state["comp_domain_input"] = "Kinase Domain"
            st.session_state["comp_lig_name_input"] = "Erlotinib"
            st.session_state["comp_smiles_input"] = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
            st.rerun()

    # Generic inputs: Blank by default unless explicitly entered in this session
    default_pocket = st.session_state.get("comp_pocket_val", "")
    default_domain = st.session_state.get("comp_domain_val", "")
    default_l_name = st.session_state.get("comp_lig_name_val", "")
    default_smiles = st.session_state.get("comp_smiles_val", "")

    # State variables
    pair_val = {"is_valid": False, "identity_pct": 0.0, "error": None}
    h_chain = h_chain_default or "A"
    m_chain = m_chain_default or "A"
    pocket_res = ""
    target_domain = ""
    final_m_res_list = []
    override_reason = ""
    ref_lig_name = ""
    ref_lig_chain = ""
    ref_lig_seq = ""

    r1_pocket_res = ""
    r2_pocket_res = ""

    # 1. Target Proteins & Binding Site Card
    with st.container(border=True):
        if workflow_mode == "matched_variant":
            st.markdown("#### 🏛️ 1️⃣ Target Receptors & Homologous Pair Validation — WT (Wild-Type Reference) vs. MT (Mutant Variant)")
            col_m1, col_m2, col_m3, col_m4 = st.columns(4)
            with col_m1:
                st.metric("WT (Wild-Type Reference)", st.session_state.get("h_id", "WT"))
                h_chain = st.selectbox(
                    "WT (Wild-Type Reference) Chain:",
                    h_chains,
                    index=h_chains.index(h_chain_default) if h_chain_default in h_chains else 0,
                    key="comp_h_chain_sel",
                )
            with col_m2:
                st.metric("MT (Mutant Variant)", st.session_state.get("m_id", "MUTANT"))
                m_chain = st.selectbox(
                    "MT (Mutant Variant) Chain:",
                    m_chains,
                    index=m_chains.index(m_chain_default) if m_chain_default in m_chains else 0,
                    key="comp_m_chain_sel",
                )
            with col_m3:
                pocket_res = st.text_input(
                    "WT (Wild-Type Reference) Pocket Residues:",
                    value=default_pocket,
                    placeholder="e.g., A:790, A:858 or B:6",
                    key="comp_pocket_input",
                )
            with col_m4:
                target_domain = st.text_input(
                    "Functional Domain:",
                    value=default_domain,
                    placeholder="e.g., Kinase Domain",
                    key="comp_domain_input",
                )
                separate_domains = st.checkbox(
                    "Domain differs between WT and MT (Strict Domain Check)",
                    value=False,
                    key="comp_domain_separate",
                    help="When enabled, each protein's domain is checked separately and ΔScore (Docking Score Difference) calculation is blocked if domains differ.",
                )
                if separate_domains:
                    mutant_target_domain = st.text_input(
                        "MT (Mutant Variant) Functional Domain:",
                        value=target_domain,
                        key="comp_mutant_domain_input",
                    )
                else:
                    mutant_target_domain = target_domain

            # Perform pair validation
            pair_val = validate_matched_pair(
                healthy_pdb_text=h_active_pdb,
                mutant_pdb_text=m_active_pdb,
                healthy_chain=h_chain,
                mutant_chain=m_chain,
                min_identity=MIN_MATCHED_PAIR_IDENTITY,
            )

            if not pair_val["is_valid"]:
                st.error(pair_val["error"])
                st.info("💡 If you are comparing two distinct proteins, switch to the 'Exploratory Cross-Protein Docking' workflow above.")
            else:
                st.success(
                    f"✅ Validated Homologous Pair (Sequence Identity: {pair_val['identity_pct']:.1f}% >= {MIN_MATCHED_PAIR_IDENTITY:.1f}%) | "
                    f"Aligned Residues: {pair_val['aligned_residue_count']} | "
                    f"Mutations: {pair_val['mutation_count']}"
                )

                # Show mapped pocket residues table
                h_res_list = [r.strip() for r in pocket_res.split(",") if r.strip()]
                if h_res_list:
                    mapped_m_res, unmapped_h_res, map_errors = map_healthy_pocket_residues(
                        h_res_list, pair_val["mapping"], m_active_pdb, m_chain, healthy_chain=h_chain
                    )

                    mapping_rows = []
                    for h_r in h_res_list:
                        m_r = pair_val["mapping"].get(h_r)
                        if not m_r and ":" not in h_r and h_chain:
                            m_r = pair_val["mapping"].get(f"{h_chain}:{h_r}")
                        status_str = "✅ Mapped" if m_r else "❌ Not found in MT (Mutant Variant)"
                        mapping_rows.append({
                            "WT (Wild-Type Reference) Residue": h_r,
                            "Mapped MT (Mutant Variant) Residue": m_r or "None",
                            "Mapping Status": status_str,
                        })

                    st.markdown("##### 🗺️ Mapped Binding-Site Residues:")
                    st.dataframe(pd.DataFrame(mapping_rows), use_container_width=True, hide_index=True)

                    if unmapped_h_res:
                        st.warning(f"⚠️ Unmapped residues not found in MT (Mutant Variant): {', '.join(unmapped_h_res)}")

                    # Manual override option with explicit reason recorded in manifest
                    override_pocket = st.checkbox(
                        "⚙️ Manual MT (Mutant Variant) Pocket Override",
                        key="chk_comp_override_pocket",
                    )
                    if override_pocket:
                        col_ov1, col_ov2 = st.columns(2)
                        with col_ov1:
                            m_pocket_override_str = st.text_input(
                                "Overridden MT (Mutant Variant) Pocket Residues:",
                                value=", ".join(mapped_m_res),
                                key="comp_override_m_pocket_input",
                            )
                        with col_ov2:
                            override_reason = st.text_input(
                                "Override Reason (required for manifest provenance):",
                                placeholder="e.g., Alternate PDB (Protein Data Bank) numbering or experimentally validated shift",
                                key="comp_override_reason_input",
                            )
                        final_m_res_list = [r.strip() for r in m_pocket_override_str.split(",") if r.strip()]
                    else:
                        final_m_res_list = mapped_m_res
                        override_reason = ""
                else:
                    final_m_res_list = []

            # Optional WT Reference Ligand expander (Section 5)
            with st.expander("🔬 Optional Protocol Validation via WT (Wild-Type Reference) Co-Crystal Ligand", expanded=False):
                st.caption("Optional: redock the co-crystallized WT reference ligand to determine empirical pose-reproduction validation status.")
                cr_col1, cr_col2, cr_col3, cr_col4 = st.columns(4)
                with cr_col1:
                    ref_lig_name = st.text_input("Reference Ligand Code (Resname):", placeholder="e.g., AQ4 or STI", key="comp_ref_lig_resname")
                with cr_col2:
                    ref_lig_chain = st.text_input("Ligand Chain:", placeholder="e.g., A", key="comp_ref_lig_chain")
                with cr_col3:
                    ref_lig_seq = st.text_input("Ligand Residue Number:", placeholder="e.g., 1001", key="comp_ref_lig_seq")
                with cr_col4:
                    ref_lig_smiles = st.text_input("Reference SMILES (optional):", placeholder="Fetched from RCSB PDB if left blank", key="comp_ref_lig_smiles")

        else:
            # Exploratory Cross-Protein Docking UI
            st.markdown("#### 🏛️ 1️⃣ Independent Receptors & Separate Pockets (Exploratory Cross-Protein)")
            st.caption("Docks any two proteins with completely independent search boxes without computing ΔScore (Docking Score Difference).")
            col_xp1, col_xp2, col_xp3, col_xp4 = st.columns(4)
            with col_xp1:
                st.metric("Receptor 1", st.session_state.get("h_id", "REC1"))
                h_chain = st.selectbox("Receptor 1 Chain:", h_chains, index=0, key="cross_r1_chain_sel")
            with col_xp2:
                st.metric("Receptor 2", st.session_state.get("m_id", "REC2"))
                m_chain = st.selectbox("Receptor 2 Chain:", m_chains, index=0, key="cross_r2_chain_sel")
            with col_xp3:
                r1_pocket_res = st.text_input(
                    "Receptor 1 Pocket Residues:",
                    value=default_pocket,
                    placeholder="e.g., A:100, A:105",
                    key="cross_r1_pocket_input",
                )
            with col_xp4:
                r2_pocket_res = st.text_input(
                    "Receptor 2 Pocket Residues:",
                    placeholder="e.g., B:50, B:55",
                    key="cross_r2_pocket_input",
                )
            target_domain = st.text_input("Description / Functional Domain:", value="Cross-protein exploratory", key="cross_domain_input")

    # 2. Comparative Docking Workflow Tabs
    tab_single_comp, tab_panel_screen = st.tabs([
        "🔬 Single Drug Comparison",
        "💊 Multi-Drug Panel Virtual Screening"
    ])

    with tab_single_comp:
        comp_svg_2d = ""
        comp_drug_props = None

        with st.container(border=True):
            st.markdown("#### 💊 2️⃣ Drug Specification & Simulation Parameters")
            # Direct PubChem search
            col_cpc1, col_cpc2 = st.columns([2.2, 1.0])
            with col_cpc1:
                comp_pubchem_query = st.text_input(
                    "🔍 Search Drug in PubChem:",
                    placeholder="e.g., Gefitinib, Osimertinib, Erlotinib",
                    key="comp_pubchem_search_input"
                )
            with col_cpc2:
                st.write("")
                st.write("")
                if st.button("🌐 Fetch Drug", key="btn_comp_pubchem_fetch", use_container_width=True):
                    if comp_pubchem_query.strip():
                        with st.spinner(f"Searching PubChem for '{comp_pubchem_query}'..."):
                            pc_res = fetch_drug_from_pubchem(comp_pubchem_query)
                            if pc_res.get("success"):
                                st.session_state["comp_lig_name_val"] = pc_res["title"]
                                st.session_state["comp_smiles_val"] = pc_res["smiles"]
                                st.session_state["comp_lig_name_input"] = pc_res["title"]
                                st.session_state["comp_smiles_input"] = pc_res["smiles"]
                                st.session_state["comp_pubchem_cid"] = {
                                    "cid": pc_res.get("cid"),
                                    "smiles": pc_res.get("smiles", ""),
                                }
                                st.success(
                                    f"✅ Fetched {pc_res['title']} (PubChem CID [Compound Identifier]: {pc_res.get('cid')})!"
                                )
                                st.rerun()
                            else:
                                st.error(pc_res.get("error", "PubChem lookup failed."))
                    else:
                        st.warning("Please enter a drug name first.")

            col_cdn1, col_cdn2 = st.columns([1, 1.6])
            with col_cdn1:
                lig_name = st.text_input("Drug Name:", value=default_l_name, key="comp_lig_name_input")
            with col_cdn2:
                lig_smiles = st.text_input("SMILES (Simplified Molecular-Input Line-Entry System):", value=default_smiles, key="comp_smiles_input")

            # 2D chemical depiction and pharmacological properties
            if lig_smiles.strip():
                tab_c_struct, tab_c_lip = st.tabs(["🧪 2D Chemical Structure", "📊 Lipinski Rules — ADME (Absorption, Distribution, Metabolism, and Excretion)"])
                with tab_c_struct:
                    comp_svg_2d = generate_2d_depiction_svg(lig_smiles, width=320, height=160)
                    if comp_svg_2d:
                        try:
                            st.image(comp_svg_2d, caption="2D (Two-Dimensional) Chemical Structure Depiction", use_container_width=True)
                        except Exception:
                            st.caption("Could not display 2D chemical depiction.")
                    else:
                        st.caption("Could not generate 2D chemical depiction for the current SMILES string.")
                with tab_c_lip:
                    comp_drug_props = calculate_drug_properties(lig_smiles)
                    if comp_drug_props and "error" not in comp_drug_props:
                        cdp1, cdp2, cdp3 = st.columns(3)
                        cdp1.metric("MW (Molecular Weight)", f"{comp_drug_props['molecular_weight']} Da (Daltons / atomic mass units)")
                        cdp2.metric("logP (Octanol–Water Partition Coefficient)", f"{comp_drug_props['logp']}")
                        cdp3.metric("TPSA (Topological Polar Surface Area)", f"{comp_drug_props['tpsa']} Å²")
                        cdp4, cdp5, cdp6 = st.columns(3)
                        cdp4.metric("HBD (Hydrogen-Bond Donors)", comp_drug_props['hbd'])
                        cdp5.metric("HBA (Hydrogen-Bond Acceptors)", comp_drug_props['hba'])
                        cdp6.metric("Rotatable Bonds", comp_drug_props['rotatable_bonds'])
                        if comp_drug_props['lipinski_pass']:
                            st.caption(
                                f"✅ Lipinski's Rule of Five Screen: {comp_drug_props.get('lipinski_violations', 0)} violation(s). "
                                "Physicochemical heuristic only, not a prediction of clinical absorption or safety."
                            )
                        else:
                            st.caption(f"⚠️ {comp_drug_props['drug_likeness']}")
                        if comp_drug_props.get("pains_alert"):
                            st.warning(
                                "⚠️ PAINS (Pan-Assay Interference Compounds) Alert: Chemotype may cause assay interference "
                                f"({comp_drug_props.get('pains_name', 'PAINS')})."
                            )
                        if comp_drug_props.get("aggregation_risk"):
                            st.warning(
                                f"⚠️ Heuristic Aggregation Risk: {comp_drug_props.get('aggregation_note', '')}"
                            )

            st.markdown("---")
            if is_advanced_mode:
                st.markdown("##### ⚙️ Advanced Comparison Settings:")
                col_mc1, col_mc2 = st.columns(2)
                with col_mc1:
                    exh = st.select_slider("Search Exhaustiveness:", options=[4, 8, 16, 32], value=8, key="comp_exh")
                with col_mc2:
                    st.caption("🔬 Identical strict protocol settings are applied to both receptors across 3 independent seeds to ensure statistical equivalence.")
            else:
                exh = 8
                st.caption("⚡ **Quick Exploratory Mode:** Search Exhaustiveness = 8 | Parallel docking across 3 independent seeds.")

        h_coordinate_source = (
            "biological_assembly_1"
            if st.session_state.get("h_assembly_pdb")
            else "deposited_or_uploaded"
        )
        m_coordinate_source = (
            "biological_assembly_1"
            if st.session_state.get("m_assembly_pdb")
            else "deposited_or_uploaded"
        )

        comp_fp_payload = {
            "workflow": workflow_mode,
            "healthy_id": st.session_state.get("h_id", ""),
            "mutant_id": st.session_state.get("m_id", ""),
            "healthy_chain": h_chain,
            "mutant_chain": m_chain,
            "healthy_coordinate_source": h_coordinate_source,
            "mutant_coordinate_source": m_coordinate_source,
            "healthy_coordinate_sha256": compute_text_sha256(h_active_pdb),
            "mutant_coordinate_sha256": compute_text_sha256(m_active_pdb),
            "pocket_res": pocket_res if workflow_mode == "matched_variant" else f"{r1_pocket_res}::{r2_pocket_res}",
            "target_domain": target_domain,
            "lig_name": lig_name,
            "lig_smiles": lig_smiles,
            "exhaustiveness": exh,
            "vina_version": vina_ver,
        }
        current_comp_fp = compute_input_fingerprint(comp_fp_payload)
        prev_comp_fp = st.session_state.get("dock_comp_input_fingerprint")
        if prev_comp_fp and prev_comp_fp != current_comp_fp:
            st.session_state.pop("matched_comparison_result", None)
            st.session_state.pop("exploratory_cross_result", None)
        st.session_state["dock_comp_input_fingerprint"] = current_comp_fp

        btn_label = (
            "🚀 Run Matched Comparison"
            if workflow_mode == "matched_variant"
            else "🚀 Run Exploratory Cross-Protein Docking"
        )

        if st.button(btn_label, type="primary", key="btn_run_comp_main"):
            if workflow_mode == "matched_variant":
                st.session_state.pop("matched_comparison_result", None)
                st.session_state.pop("exploratory_cross_result", None)

                h_res_list = [r.strip() for r in pocket_res.split(",") if r.strip()]
                if not h_res_list:
                    st.error("Docking was not run: define a validated binding site.")
                    return
                if not target_domain.strip():
                    st.error("Please enter the required functional domain.")
                    return
                if not lig_name.strip():
                    st.error("Please enter the drug name.")
                    return
                if not lig_smiles.strip():
                    st.error("Please enter the drug SMILES (Simplified Molecular-Input Line-Entry System) string.")
                    return

                # Strict pair validation check
                if not pair_val.get("is_valid"):
                    st.error(pair_val.get("error") or "Selected chains are not a validated WT/variant pair.")
                    return

                # Check manual override requirements
                if st.session_state.get("chk_comp_override_pocket") and not override_reason.strip():
                    st.error("Please provide a reason for the manual MT (Mutant Variant) pocket override for manifest provenance.")
                    return

                if not final_m_res_list:
                    st.error("MT (Mutant Variant) binding pocket is undefined or no homologous residues were mapped.")
                    return

                # Optional WT reference ligand extraction
                wt_ref_meta = None
                if ref_lig_name.strip():
                    try:
                        _, ref_pdb_text = extract_reference_ligand(
                            h_active_pdb,
                            ref_lig_name.strip(),
                            chain_id=ref_lig_chain.strip() or None,
                            resseq=ref_lig_seq.strip() or None,
                        )
                        wt_ref_meta = {
                            "resname": ref_lig_name.strip(),
                            "chain": ref_lig_chain.strip(),
                            "resseq": ref_lig_seq.strip(),
                            "ref_ligand_pdb_text": ref_pdb_text,
                            "ref_smiles": ref_lig_smiles.strip() if ref_lig_smiles.strip() else None,
                        }
                    except Exception as ex_ref:
                        st.warning(f"Could not extract reference ligand: {ex_ref}")

                pair_validation_meta = {
                    "is_valid": pair_val["is_valid"],
                    "identity_pct": pair_val["identity_pct"],
                    "min_identity": MIN_MATCHED_PAIR_IDENTITY,
                    "healthy_chain": h_chain,
                    "mutant_chain": m_chain,
                    "aligned_residue_count": pair_val["aligned_residue_count"],
                    "mutation_count": pair_val["mutation_count"],
                    "mutations": pair_val["mutations"],
                    "mapping": pair_val["mapping"],
                    "healthy_pocket_residues": h_res_list,
                    "mutant_pocket_residues": final_m_res_list,
                    "override_reason": override_reason,
                }

                with st.spinner("Running matched WT (Wild-Type Reference) vs. MT (Mutant Variant) docking across 3 independent seeds..."):
                    h_meta = {"pdb_id": st.session_state.get("h_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
                    m_meta = {"pdb_id": st.session_state.get("m_id"), "target_domain": mutant_target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}

                    try:
                        comp_result = run_matched_docking_comparison(
                            healthy_pdb_text=h_active_pdb,
                            mutant_pdb_text=m_active_pdb,
                            healthy_target_meta=h_meta,
                            mutant_target_meta=m_meta,
                            ligand_smiles=lig_smiles,
                            ligand_name=lig_name,
                            healthy_pocket_residues=h_res_list,
                            mutant_pocket_residues=final_m_res_list,
                            target_domain=target_domain,
                            exhaustiveness=exh,
                            input_fingerprint=current_comp_fp,
                            pair_validation_meta=pair_validation_meta,
                            wt_reference_ligand_meta=wt_ref_meta,
                            pubchem_cid=_get_pubchem_cid(lig_smiles, "comp_pubchem_cid"),
                            keep_waters=keep_crystal_waters,
                        )
                        comp_result["input_fingerprint"] = current_comp_fp
                        st.session_state["matched_comparison_result"] = comp_result
                        st.success("✅ Matched comparison completed successfully!")
                    except Exception as e:
                        st.error(f"❌ Matched comparison failed: {e}")

            else:
                # Exploratory cross-protein docking
                st.session_state.pop("matched_comparison_result", None)
                st.session_state.pop("exploratory_cross_result", None)

                r1_res_list = [r.strip() for r in r1_pocket_res.split(",") if r.strip()]
                r2_res_list = [r.strip() for r in r2_pocket_res.split(",") if r.strip()]

                if not r1_res_list or not r2_res_list:
                    st.error("Exploratory docking requires pocket residues defined for both receptors.")
                    return
                if not lig_name.strip():
                    st.error("Please enter the drug name.")
                    return
                if not lig_smiles.strip():
                    st.error("Please enter the drug SMILES (Simplified Molecular-Input Line-Entry System) string.")
                    return

                # Validate pocket chains match selected chains (Section 4)
                try:
                    validate_pocket_chain(r1_res_list, h_chain, "Receptor 1")
                    validate_pocket_chain(r2_res_list, m_chain, "Receptor 2")
                except ValueError as chain_err:
                    st.error(f"❌ Selected chain mismatch: {chain_err}")
                    return

                with st.spinner("Running independent exploratory cross-protein docking across 3 independent seeds..."):
                    r1_meta = {
                        "pdb_id": st.session_state.get("h_id", "REC1"),
                        "selected_chain": h_chain,
                        "target_domain": target_domain,
                        "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    }
                    r2_meta = {
                        "pdb_id": st.session_state.get("m_id", "REC2"),
                        "selected_chain": m_chain,
                        "target_domain": target_domain,
                        "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    }

                    try:
                        cross_result = run_exploratory_cross_docking(
                            receptor1_pdb_text=h_active_pdb,
                            receptor2_pdb_text=m_active_pdb,
                            receptor1_target_meta=r1_meta,
                            receptor2_target_meta=r2_meta,
                            ligand_smiles=lig_smiles,
                            ligand_name=lig_name,
                            receptor1_pocket_residues=r1_res_list,
                            receptor2_pocket_residues=r2_res_list,
                            exhaustiveness=exh,
                            input_fingerprint=current_comp_fp,
                            pubchem_cid=_get_pubchem_cid(lig_smiles, "comp_pubchem_cid"),
                            keep_waters=keep_crystal_waters,
                        )
                        st.session_state["exploratory_cross_result"] = cross_result
                        st.success("✅ Exploratory cross-protein docking completed successfully!")
                    except Exception as e:
                        st.error(f"❌ Exploratory cross-protein docking failed: {e}")

        # Render Matched Results
        comp = st.session_state.get("matched_comparison_result")
        if workflow_mode == "matched_variant" and comp:
            st.markdown("---")
            st.subheader(f"📊 Comparative Docking Results: {comp['comparison_label']}")

            # Validation status badge and reason
            val_status = comp.get("validation_status", "exploratory—no reference validation")
            val_reason = comp.get("validation_reason", "")
            val_norm = str(val_status).strip().lower()
            if val_norm == "protocol validated on wt reference":
                st.success("🏆 **Scientific Validation Status:** Protocol validated on WT (Wild-Type Reference) reference")
            elif val_norm in ("exploratory—no reference validation", "exploratory: no pose validation"):
                st.info("🧪 **Scientific Validation Status:** Exploratory: no reference validation")
            else:
                st.warning("⚠️ **Scientific Validation Status:** Not validated for affinity comparison")

            if val_reason:
                st.caption(f"📝 **Validation Rationale:** {val_reason}")
            mutant_val_status = comp.get("mutant_validation_status")
            if mutant_val_status and mutant_val_status != val_status:
                st.caption(
                    f"🧬 **MT (Mutant Variant) Arm Status:** {mutant_val_status} — "
                    f"{comp.get('mutant_validation_reason', '')}"
                )

            redocking_val = comp.get("redocking_validation")
            if redocking_val and "top1" in redocking_val and "recovery" in redocking_val:
                t1 = redocking_val["top1"]
                rec = redocking_val["recovery"]
                cutoff = redocking_val.get("rmsd_cutoff", 2.0)
                with st.expander("🔬 WT (Wild-Type Reference) Redocking Validation Metrics — RMSD"):
                    st.markdown(
                        f"- **Acceptance Cutoff — RMSD:** `<= {cutoff} Å`\n"
                        f"- **Top-1 Vina pose RMSD:** `{t1['rmsd_angstrom']:.3f} Å` — Score: `{t1['score']} kcal/mol`\n"
                        f"- **Best Top-10 recovery RMSD:** `{rec['rmsd_angstrom']:.3f} Å` — Seed `{rec['seed']}`, Rank `{rec['rank']}`, Score: `{rec['score']} kcal/mol`\n"
                        f"- **Independent seed support for cluster:** `{rec['cluster_seed_count']}` seeds `{rec['supporting_seeds']}`\n"
                        f"- **Mapping method:** `{redocking_val.get('method', 'symmetry-aware heavy-atom mapped RMSD (no-fit)')}`"
                    )

            # Protonation assumption note (Section 6)
            st.info(f"ℹ️ {comp.get('protonation_note', PROTONATION_ASSUMPTION_NOTE)}")

            # Calculate pharmacological interactions
            h_clean_text = ""
            m_clean_text = ""
            try:
                with open(comp["healthy"]["clean_pdb"], "r", encoding="utf-8") as f:
                    h_clean_text = f.read()
                with open(comp["mutant"]["clean_pdb"], "r", encoding="utf-8") as f:
                    m_clean_text = f.read()
            except Exception:
                pass

            # Use consensus dominant cluster pose (Section 4)
            h_dom_pose = comp["healthy"].get("dominant_cluster_pose")
            m_dom_pose = comp["mutant"].get("dominant_cluster_pose")
            h_top_pose = h_dom_pose.pdbqt_block if h_dom_pose else ""
            m_top_pose = m_dom_pose.pdbqt_block if m_dom_pose else ""

            h_inter = analyze_protein_ligand_interactions(h_clean_text, h_top_pose)
            m_inter = analyze_protein_ligand_interactions(m_clean_text, m_top_pose)

            h_pocket_list = comp["healthy"].get("pocket_residues", [])
            m_pocket_list = comp["mutant"].get("pocket_residues", [])

            comp_inter = compare_interactions(h_inter, m_inter, mutation_residue_nums=h_pocket_list)
            delta_score = comp.get("delta_score")
            score_diff_info = (
                interpret_score_difference(
                    delta_score,
                    comp.get("delta_score_uncertainty_kcal"),
                )
                if delta_score is not None
                else None
            )

            # Top KPI Metric Cards
            with st.container(border=True):
                kpi_c1, kpi_c2, kpi_c3, kpi_c4 = st.columns(4)
                kpi_c1.metric("WT (Wild-Type Reference) Median Score", f"{comp['healthy']['dominant_cluster_median']} kcal/mol")
                kpi_c2.metric("MT (Mutant Variant) Median Score", f"{comp['mutant']['dominant_cluster_median']} kcal/mol")
                kpi_c3.metric(
                    "ΔScore",
                    (
                        f"{delta_score:+.2f} kcal/mol"
                        if delta_score is not None
                        else "Not evaluated"
                    ),
                )
                kpi_c4.metric(
                    "Score-Difference Interpretation",
                    score_diff_info["category"] if score_diff_info else "Not evaluated",
                )

            col_s1, col_s2 = st.columns(2)
            with col_s1:
                st.write(
                    f"**WT (Wild-Type Reference) Supporting Seeds:** {comp['healthy']['dominant_cluster_seeds']} "
                    f"({comp['healthy']['seed_count']}) | Consensus Cluster: #{comp['healthy']['dominant_cluster_id']}"
                )
            with col_s2:
                st.write(
                    f"**MT (Mutant Variant) Supporting Seeds:** {comp['mutant']['dominant_cluster_seeds']} "
                    f"({comp['mutant']['seed_count']}) | Consensus Cluster: #{comp['mutant']['dominant_cluster_id']}"
                )

            if delta_score is not None:
                st.caption(
                    "Note: ΔScore (Docking Score Difference) = MT (Mutant Variant) Median Score − WT (Wild-Type Reference) Median Score "
                    "(a more negative value means a lower — more favorable — Vina score under this protocol)."
                )
                delta_uncertainty = comp.get("delta_score_uncertainty_kcal")
                if delta_uncertainty is not None:
                    st.caption(
                        f"Estimated ΔScore uncertainty from per-seed best-score variability: "
                        f"±{delta_uncertainty:.2f} kcal/mol. Differences smaller than "
                        f"max(±{delta_uncertainty:.2f}, 0.50) kcal/mol are reported as within "
                        "protocol noise, with no direction claimed."
                    )
            else:
                st.warning(
                    "ΔScore (Docking Score Difference) was not evaluated: the binding-mode "
                    "gate requires the WT and MT consensus poses to represent the same site "
                    "and geometry. The per-arm median scores above are shown for transparency "
                    "only and are not a matched comparison result."
                )

            # Organized comparison result tabs
            comp_tab_3d, comp_tab_resist = st.tabs([
                "🔮 Synchronized 3D (Three-Dimensional) Comparison — WT vs. MT",
                "🔬 Exploratory Contact Candidates"
            ])

            with comp_tab_3d:
                with st.container(border=True):
                    c_col1, c_col2, c_col3 = st.columns(3)
                    with c_col1:
                        comp_show_ligand = st.checkbox("Show Ligand", value=True, key="dock_comp_show_ligand")
                        if comp_show_ligand:
                            comp_lig_style = st.selectbox("Ligand Style:", ["Sticks", "Ball & Stick", "Spheres"], key="dock_comp_lig_style")
                        else:
                            comp_lig_style = "Sticks"
                    with c_col2:
                        comp_show_grid = st.checkbox("📦 Search Grid Box", value=False, key="dock_comp_show_grid")
                        comp_show_surface = st.checkbox("Show Surface", value=False, key="dock_comp_show_surface")
                    with c_col3:
                        if comp_show_surface:
                            comp_surf_type = st.selectbox("Surface Type:", ["MS (Molecular Surface)", "SAS (Solvent-Accessible Surface)", "VDW (van der Waals Surface)"], key="dock_comp_surf_type")
                            comp_surf_opacity = st.slider("Surface Opacity:", 0.0, 1.0, 0.3, key="dock_comp_surf_opacity")
                        else:
                            comp_surf_type = "MS"
                            comp_surf_opacity = 0.3

                v_col1, v_col2 = st.columns(2)
                with v_col1:
                    st.markdown(
                        f"**🟢 WT (Wild-Type Reference): {comp['healthy']['target_id']}** "
                        f"[Cluster #{comp['healthy']['dominant_cluster_id']} | "
                        f"Seeds: {comp['healthy']['seed_count']} | "
                        f"Median: {comp['healthy']['dominant_cluster_median']} kcal/mol]"
                    )
                    st.caption("🟣 Magenta: WT (Wild-Type Reference) Polar-Contact Residues | 🟠 Orange: Binding-Pocket Residues")
                    h_view_html = render_docked_complex_3d(
                        comp['healthy']['clean_pdb'],
                        h_top_pose,
                        "",
                        show_surface=comp_show_surface,
                        surface_opacity=comp_surf_opacity,
                        surface_type=comp_surf_type.split()[0],
                        pocket_residues=h_pocket_list,
                        hbonds=h_inter.get("polar_contact_candidates") or h_inter.get("hbonds"),
                        show_ligand=comp_show_ligand,
                        ligand_style=comp_lig_style,
                        grid_box=comp['healthy'].get('grid'),
                        show_grid_box=comp_show_grid
                    )
                    st.iframe(h_view_html, height=480)

                with v_col2:
                    st.markdown(
                        f"**🔴 MT (Mutant Variant): {comp['mutant']['target_id']}** "
                        f"[Cluster #{comp['mutant']['dominant_cluster_id']} | "
                        f"Seeds: {comp['mutant']['seed_count']} | "
                        f"Median: {comp['mutant']['dominant_cluster_median']} kcal/mol]"
                    )
                    st.caption("🟣 Magenta: MT (Mutant Variant) Polar-Contact Residues | 🟠 Orange: Binding-Pocket Residues")
                    m_view_html = render_docked_complex_3d(
                        comp['mutant']['clean_pdb'],
                        m_top_pose,
                        "",
                        show_surface=comp_show_surface,
                        surface_opacity=comp_surf_opacity,
                        surface_type=comp_surf_type.split()[0],
                        pocket_residues=m_pocket_list,
                        hbonds=m_inter.get("polar_contact_candidates") or m_inter.get("hbonds"),
                        show_ligand=comp_show_ligand,
                        ligand_style=comp_lig_style,
                        grid_box=comp['mutant'].get('grid'),
                        show_grid_box=comp_show_grid
                    )
                    st.iframe(m_view_html, height=480)

            with comp_tab_resist:
                if score_diff_info is None:
                    st.warning("### ΔScore Not Evaluated (Binding-Mode Gate)")
                    st.write(
                        "The WT and MT consensus poses do not share the same binding mode, "
                        "so a docking-score difference would be uninterpretable. No direction "
                        "or resistance claim is made."
                    )
                elif score_diff_info["direction"] == "less favorable":
                    st.warning(f"### {score_diff_info['category']}")
                    st.write(score_diff_info["summary"])
                elif score_diff_info["direction"] == "more favorable":
                    st.success(f"### {score_diff_info['category']}")
                    st.write(score_diff_info["summary"])
                else:
                    st.info(f"### {score_diff_info['category']}")
                    st.write(score_diff_info["summary"])
                st.info(CONTACTS_DISCLAIMER)

                with st.container(border=True):
                    rc_col1, rc_col2, rc_col3 = st.columns(3)
                    p_healthy = comp_inter.get("healthy_polar_contact_count", comp_inter.get("healthy_polar_count", 0))
                    p_mutant = comp_inter.get("mutant_polar_contact_count", comp_inter.get("mutant_polar_count", 0))
                    lost_polar = comp_inter.get("lost_polar_contact_residues", comp_inter.get("lost_polar_residues", []))
                    gained_polar = comp_inter.get("gained_polar_contact_residues", comp_inter.get("gained_polar_residues", []))
                    conserved_polar = comp_inter.get("conserved_polar_contact_residues", comp_inter.get("conserved_polar_residues", []))

                    rc_col1.metric("WT (Wild-Type Reference) Polar Contacts", p_healthy)
                    rc_col2.metric("MT (Mutant Variant) Polar Contacts", p_mutant)
                    rc_col3.metric("Lost Polar Contacts in MT", len(lost_polar))

                    if lost_polar:
                        st.warning(f"⚠️ **Polar contact candidates lost in MT (Mutant Variant):** {', '.join(lost_polar)}")
                    if gained_polar:
                        st.success(f"✨ **Polar contact candidates gained in MT (Mutant Variant):** {', '.join(gained_polar)}")
                    if conserved_polar:
                        st.info(f"🛡️ **Conserved polar contact candidates:** {', '.join(conserved_polar)}")

                    col_tb1, col_tb2 = st.columns(2)
                    with col_tb1:
                        st.markdown("**🟢 WT (Wild-Type Reference) Polar Contact Candidates:**")
                        h_polar = h_inter.get("polar_contact_candidates") or h_inter.get("hbonds", [])
                        if h_polar:
                            df_h_hb = pd.DataFrame(h_polar).rename(columns={
                                "residue": "Residue", "chain": "Chain", "rec_atom": "Receptor Atom", "lig_atom": "Ligand Atom", "distance_angstrom": "Distance — Å"
                            })[["Residue", "Chain", "Receptor Atom", "Ligand Atom", "Distance — Å"]]
                            st.dataframe(df_h_hb, use_container_width=True, hide_index=True)
                        else:
                            st.caption("No polar contact candidates within <= 3.5 Å.")

                    with col_tb2:
                        st.markdown("**🔴 MT (Mutant Variant) Polar Contact Candidates:**")
                        m_polar = m_inter.get("polar_contact_candidates") or m_inter.get("hbonds", [])
                        if m_polar:
                            df_m_hb = pd.DataFrame(m_polar).rename(columns={
                                "residue": "Residue", "chain": "Chain", "rec_atom": "Receptor Atom", "lig_atom": "Ligand Atom", "distance_angstrom": "Distance — Å"
                            })[["Residue", "Chain", "Receptor Atom", "Ligand Atom", "Distance — Å"]]
                            st.dataframe(df_m_hb, use_container_width=True, hide_index=True)
                        else:
                            st.caption("No polar contact candidates within <= 3.5 Å.")

        # Render Exploratory Cross Results
        cross_res = st.session_state.get("exploratory_cross_result")
        if workflow_mode == "exploratory_cross" and cross_res:
            st.markdown("---")
            st.subheader(f"📊 {cross_res['comparison_label']}")
            st.info(f"ℹ️ {cross_res.get('protonation_note', PROTONATION_ASSUMPTION_NOTE)}")
            st.warning(
                "⚠️ Strict Scientific Notice: ΔScore (Docking Score Difference) and comparative ranking "
                "are not computed across non-homologous proteins or unmatched pockets."
            )

            col_xp_res1, col_xp_res2 = st.columns(2)
            with col_xp_res1:
                r1_data = cross_res["receptor1"]
                with st.container(border=True):
                    st.markdown(f"### 🏛️ Receptor 1: {r1_data['target_id']}")
                    if r1_data.get("cluster_reproducible"):
                        st.success("✅ Reproducible exploratory cluster (>= 2 independent seeds)")
                        st.metric("Median Score", f"{r1_data['dominant_cluster_median']} kcal/mol")
                    else:
                        st.warning("⚠️ Exploratory pose without independent-seed support")
                        r1_dom_p = r1_data.get("dominant_cluster_pose")
                        r1_disp_score = r1_dom_p.score if r1_dom_p else r1_data.get("dominant_cluster_median")
                        st.metric("Exploratory Pose Score", f"{r1_disp_score} kcal/mol")

                    r1_dom_pose = r1_data.get("dominant_cluster_pose")
                    if r1_dom_pose:
                        st.caption(
                            f"**Cluster:** #{r1_data['dominant_cluster_id']} | "
                            f"**Pose:** Rank {r1_dom_pose.rank} (Seed {r1_dom_pose.seed}, Score {r1_dom_pose.score} kcal/mol) | "
                            f"**Supporting Seeds:** {r1_data['dominant_cluster_seeds']}"
                        )
                    else:
                        st.caption(f"**Cluster:** #{r1_data['dominant_cluster_id']} | **Supporting Seeds:** {r1_data['dominant_cluster_seeds']}")

                    r1_pose = r1_dom_pose.pdbqt_block if r1_dom_pose else ""
                    r1_html = render_docked_complex_3d(
                        r1_data["clean_pdb"],
                        r1_pose,
                        "",
                        pocket_residues=r1_data["pocket_residues"],
                        grid_box=r1_data.get("grid"),
                    )
                    st.iframe(r1_html, height=400)

            with col_xp_res2:
                r2_data = cross_res["receptor2"]
                with st.container(border=True):
                    st.markdown(f"### 🏛️ Receptor 2: {r2_data['target_id']}")
                    if r2_data.get("cluster_reproducible"):
                        st.success("✅ Reproducible exploratory cluster (>= 2 independent seeds)")
                        st.metric("Median Score", f"{r2_data['dominant_cluster_median']} kcal/mol")
                    else:
                        st.warning("⚠️ Exploratory pose without independent-seed support")
                        r2_dom_p = r2_data.get("dominant_cluster_pose")
                        r2_disp_score = r2_dom_p.score if r2_dom_p else r2_data.get("dominant_cluster_median")
                        st.metric("Exploratory Pose Score", f"{r2_disp_score} kcal/mol")

                    r2_dom_pose = r2_data.get("dominant_cluster_pose")
                    if r2_dom_pose:
                        st.caption(
                            f"**Cluster:** #{r2_data['dominant_cluster_id']} | "
                            f"**Pose:** Rank {r2_dom_pose.rank} (Seed {r2_dom_pose.seed}, Score {r2_dom_pose.score} kcal/mol) | "
                            f"**Supporting Seeds:** {r2_data['dominant_cluster_seeds']}"
                        )
                    else:
                        st.caption(f"**Cluster:** #{r2_data['dominant_cluster_id']} | **Supporting Seeds:** {r2_data['dominant_cluster_seeds']}")

                    r2_pose = r2_dom_pose.pdbqt_block if r2_dom_pose else ""
                    r2_html = render_docked_complex_3d(
                        r2_data["clean_pdb"],
                        r2_pose,
                        "",
                        pocket_residues=r2_data["pocket_residues"],
                        grid_box=r2_data.get("grid"),
                    )
                    st.iframe(r2_html, height=400)

    # 3. Multi-Drug Panel Virtual Screening Tab
    with tab_panel_screen:
        with st.container(border=True):
            st.markdown("#### 💊 Multi-Drug Virtual Screening & Comparison")
            st.markdown("Compare multiple compounds against the mutation in a single exploratory screening batch to compute docking score differences:")
            screening_panels = {
                "EGFR (Epidermal Growth Factor Receptor) TKI (Tyrosine Kinase Inhibitor) Panel — NSCLC (Non-Small Cell Lung Cancer)": [
                    ("Erlotinib (1st Gen)", "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"),
                    ("Gefitinib (1st Gen)", "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4"),
                    ("Osimertinib (3rd Gen)", "CN1CCN(CC1)C2=CC(=C(C=C2)NC3=NC=CC(=N3)C4=CN(C5=CC=CC=C54)C)NC(=O)C=C")
                ],
                "BCR-ABL Inhibitors — CML (Chronic Myeloid Leukemia)": [
                    ("Imatinib (1st Gen)", "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"),
                    ("Dasatinib (2nd Gen)", "CC1=C(C(=CC=C1)Cl)NC(=O)C2=CN=C(S2)NC3=CC(=NC(=N3)C)N4CCN(CC4)CCO")
                ],
                "Clinical Antivirals": [
                    ("Nirmatrelvir (Paxlovid)", "CC(C)(C)C(C(=O)NC(C1CC1)C#N)NC(=O)C2CC3(C2)C(C(O3)(F)F)C(=O)O"),
                    ("Remdesivir", "CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)(C#N)C2=CC=C3N2N=CN=C3N)O)O)OC4=CC=CC=C4")
                ]
            }
            chosen_panel_name = st.selectbox("Select Screening Drug Panel:", list(screening_panels.keys()), key="screen_panel_select")
            chosen_drugs = screening_panels[chosen_panel_name]

            if workflow_mode != "matched_variant":
                st.info("💡 Panel virtual screening in this tab is designed for matched WT (Wild-Type Reference) vs. MT (Mutant Variant) pairs. For exploratory cross-protein docking, screen drugs individually in the first tab.")
            else:
                screen_fp_payload = {
                    "workflow": "screening_panel",
                    "healthy_id": st.session_state.get("h_id", ""),
                    "mutant_id": st.session_state.get("m_id", ""),
                    "healthy_chain": h_chain,
                    "mutant_chain": m_chain,
                    "healthy_coordinate_source": h_coordinate_source,
                    "mutant_coordinate_source": m_coordinate_source,
                    "healthy_coordinate_sha256": compute_text_sha256(h_active_pdb),
                    "mutant_coordinate_sha256": compute_text_sha256(m_active_pdb),
                    "pocket_res": pocket_res,
                    "target_domain": target_domain,
                    "chosen_panel_name": chosen_panel_name,
                    "vina_version": vina_ver,
                }
                current_screen_fp = compute_input_fingerprint(screen_fp_payload)
                prev_screen_fp = st.session_state.get("dock_screen_input_fingerprint")
                if prev_screen_fp and prev_screen_fp != current_screen_fp:
                    st.session_state.pop("screening_results", None)
                st.session_state["dock_screen_input_fingerprint"] = current_screen_fp

                if st.button("🚀 Run Panel Virtual Screening", key="btn_run_panel_screen", type="primary"):
                    st.session_state.pop("screening_results", None)
                    screen_pocket_list = [r.strip() for r in pocket_res.split(",") if r.strip()]
                    if not screen_pocket_list:
                        st.error("Please define binding-pocket residues in the upper panel first.")
                    elif not pair_val.get("is_valid"):
                        st.error(pair_val.get("error") or "Selected chains are not a validated WT/variant pair.")
                    elif not final_m_res_list:
                        st.error("MT (Mutant Variant) binding pocket is undefined or no homologous residues were mapped.")
                    else:
                        progress_bar = st.progress(0.0)
                        screen_results = []
                        h_meta_s = {"pdb_id": st.session_state.get("h_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
                        m_meta_s = {"pdb_id": st.session_state.get("m_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}

                        for idx, (d_name, d_smiles) in enumerate(chosen_drugs):
                            try:
                                drug_fp_payload = {
                                    "workflow": "screening_panel_drug",
                                    "healthy_id": st.session_state.get("h_id", ""),
                                    "mutant_id": st.session_state.get("m_id", ""),
                                    "healthy_chain": h_chain,
                                    "mutant_chain": m_chain,
                                    "healthy_coordinate_source": h_coordinate_source,
                                    "mutant_coordinate_source": m_coordinate_source,
                                    "healthy_coordinate_sha256": compute_text_sha256(h_active_pdb),
                                    "mutant_coordinate_sha256": compute_text_sha256(m_active_pdb),
                                    "pocket_res": pocket_res,
                                    "target_domain": target_domain,
                                    "vina_version": vina_ver,
                                    "exhaustiveness": 4,
                                    "drug_name": d_name,
                                    "drug_smiles": d_smiles,
                                }
                                drug_input_fp = compute_input_fingerprint(drug_fp_payload)
                                cr = run_matched_docking_comparison(
                                    healthy_pdb_text=h_active_pdb,
                                    mutant_pdb_text=m_active_pdb,
                                    healthy_target_meta=h_meta_s,
                                    mutant_target_meta=m_meta_s,
                                    ligand_smiles=d_smiles,
                                    ligand_name=d_name,
                                    healthy_pocket_residues=screen_pocket_list,
                                    mutant_pocket_residues=final_m_res_list,
                                    target_domain=target_domain,
                                    exhaustiveness=4,
                                    input_fingerprint=drug_input_fp,
                                    pair_validation_meta=(
                                        {
                                            "is_valid": bool(pair_val.get("is_valid")),
                                            "identity_pct": pair_val.get("identity_pct"),
                                        }
                                        if pair_val.get("is_valid")
                                        else None
                                    ),
                                    keep_waters=keep_crystal_waters,
                                )
                                wt_s = cr['healthy']['dominant_cluster_median']
                                mt_s = cr['mutant']['dominant_cluster_median']
                                delta_s = cr.get('delta_score')
                                noise = cr.get('delta_score_uncertainty_kcal')
                                neutral_band = (
                                    max(float(noise), 0.5)
                                    if noise is not None
                                    else SCORE_DIFF_NEUTRAL_BAND_KCAL
                                )
                                moderate_band = max(SCORE_DIFF_MODERATE_BAND_KCAL, neutral_band)
                                if delta_s is None:
                                    r_tag = "Not evaluated (binding-mode gate)"
                                elif delta_s >= moderate_band:
                                    r_tag = f"Score increase >= +{moderate_band:.2f} kcal/mol"
                                elif delta_s >= neutral_band:
                                    r_tag = (
                                        f"Score increase +{neutral_band:.2f} to +{moderate_band:.2f} kcal/mol"
                                    )
                                elif delta_s <= -moderate_band:
                                    r_tag = f"Score decrease <= -{moderate_band:.2f} kcal/mol"
                                else:
                                    r_tag = (
                                        f"Within protocol noise — |ΔScore| < {neutral_band:.2f} kcal/mol"
                                    )

                                screen_results.append({
                                    "Drug": d_name,
                                    "WT (Wild-Type Reference) kcal/mol": wt_s,
                                    "MT (Mutant Variant) kcal/mol": mt_s,
                                    "ΔScore": delta_s,
                                    "Exploratory Score-Difference Flag": r_tag
                                })
                            except Exception as ex:
                                st.warning(f"Could not screen {d_name}: {ex}")
                            progress_bar.progress((idx + 1) / len(chosen_drugs))

                        if screen_results:
                            st.session_state["screening_results"] = screen_results
                            st.success("✅ Virtual panel screening completed successfully!")

                saved_screen_res = st.session_state.get("screening_results")
                if saved_screen_res:
                    df_screen = pd.DataFrame(saved_screen_res)
                    st.dataframe(df_screen, use_container_width=True, hide_index=True)

                    fig_screen = go.Figure()
                    d_names = [r["Drug"] for r in saved_screen_res]
                    wt_vals = [r["WT (Wild-Type Reference) kcal/mol"] for r in saved_screen_res]
                    mt_vals = [r["MT (Mutant Variant) kcal/mol"] for r in saved_screen_res]

                    fig_screen.add_trace(go.Bar(
                        x=d_names, y=wt_vals, name="WT (Wild-Type Reference)",
                        marker_color="#2ea043"
                    ))
                    fig_screen.add_trace(go.Bar(
                        x=d_names, y=mt_vals, name="MT (Mutant Variant)",
                        marker_color="#da3633"
                    ))
                    fig_screen.update_layout(
                        title="Computational Docking Score Comparison Across Drug Panel — kcal/mol (Lower Is Better)",
                        xaxis_title="Drug Compound",
                        yaxis_title="Docking Score — kcal/mol",
                        barmode="group",
                        template="plotly_dark",
                        height=380,
                        margin=dict(l=20, r=20, t=40, b=20)
                    )
                    st.plotly_chart(fig_screen, use_container_width=True)

