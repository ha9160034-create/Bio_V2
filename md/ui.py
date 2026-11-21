# Molecular Dynamics Streamlit UI and analysis dashboard

from __future__ import annotations

import re
import json
import shutil
import zipfile
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, List

import streamlit as st
import streamlit.components.v1 as components
import py3Dmol

from md.engine import detect_gromacs, get_gmx_version, check_backend_tools, generate_standalone_scripts
from md.topology import (
    TopologyBuildError,
    extract_best_pose_from_docking_run,
    assemble_protein_ligand_complex,
    prepare_complex_topology,
    sanitize_resname,
)
from md.simulation import (
    SimulationConfig,
    SimulationError,
    BiophysicalMDSimulator,
    create_md_run_directory,
    compute_rhombic_dodecahedron_volume,
    generate_mdp_ions,
    generate_mdp_minimization,
    generate_mdp_nve,
    generate_mdp_nvt,
    generate_mdp_npt,
    generate_mdp_production,
    parse_hbond_num_xvg,
)
from md.analysis import (
    RMSDRecord,
    RMSFRecord,
    MMPBSA_NOT_AVAILABLE_MESSAGE,
    parse_xvg_file,
    load_mmpbsa_results,
    summarize_timeseries,
    load_last_frame_pdb,
    create_rmsd_plotly,
    create_rmsf_plotly,
    create_hbond_count_plotly,
    create_mmpbsa_breakdown_plotly,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DOCKING_RUNS_DIR = PROJECT_ROOT / "docking_runs"
MD_RUNS_DIR = PROJECT_ROOT / "md_runs"

DISCLAIMER_MD = (
    "MD and MM-PB(GB)SA results are computational predictions for relative ranking, "
    "not experimental measurements."
)


def render_mol_3d(
    pdb_text: str,
    height: int = 400,
    ligand_resname: str = "LIG",
    show_surface: bool = True,
    surface_opacity: float = 0.5,
    surface_type: str = "MS",
    surface_color: str = "#ECEFF1",
):
    """Render a 3D interactive molecular visualization using py3Dmol."""
    viewer = py3Dmol.view(width="100%", height=height)
    viewer.addModel(pdb_text, "pdb")
    viewer.setStyle({"hetflag": False}, {"cartoon": {"color": "spectrum", "opacity": 0.85}})

    surf_code = surface_type.split()[0] if surface_type else "MS"
    surf_kind = (
        py3Dmol.MS
        if (surface_type == "MS" or surf_code == "MS")
        else (py3Dmol.VDW if (surface_type == "VDW" or surf_code == "VDW") else py3Dmol.SAS)
    )
    if show_surface:
        viewer.addSurface(
            surf_kind,
            {"opacity": surface_opacity, "color": surface_color},
            {"hetflag": False},
        )

    ligand_style = {
        "stick": {"colorscheme": "greenCarbon", "radius": 0.28},
        "sphere": {"colorscheme": "greenCarbon", "radius": 0.5},
    }
    viewer.setStyle({"hetflag": True}, ligand_style)
    viewer.setStyle(
        {"resn": [ligand_resname, "JZ4", "UNL"]},
        ligand_style,
    )
    viewer.zoomTo()

    style_tag = (
        "<style>html, body { margin: 0 !important; padding: 0 !important; "
        "overflow: hidden !important; width: 100% !important; height: 100% !important; } "
        "div[id^='3dmolviewer'] { overflow: hidden !important; width: 100% !important; height: 100% !important; }</style>"
    )
    viewer.startjs = style_tag + viewer.startjs
    components.html(viewer._make_html(), height=height, scrolling=False)


def scan_available_docking_runs() -> List[Dict[str, Any]]:
    """Scan docking_runs/ for completed runs with manifests and poses."""
    runs: List[Dict[str, Any]] = []
    if not DOCKING_RUNS_DIR.exists():
        return runs

    for d in sorted(DOCKING_RUNS_DIR.iterdir(), reverse=True):
        if not d.is_dir() or d.name.startswith("."):
            continue
        has_poses = bool(list(d.glob("out_seed_*.pdbqt")) or list(d.glob("*_out.pdbqt")))
        has_rec = (d / "receptor_clean.pdb").is_file() or bool(list(d.glob("*.pdb")))
        if not (has_poses and has_rec):
            continue

        target_id = "Unknown"
        lig_name = "Ligand"
        top_score = None
        val_stat = "exploratory"
        manifest_file = d / "manifest.json"
        if manifest_file.is_file():
            try:
                m = json.loads(manifest_file.read_text(encoding="utf-8"))
                target_id = m.get("target", {}).get("pdb_id", target_id)
                lig_name = m.get("ligand", {}).get("ligand_name") or m.get("ligand", {}).get("name", lig_name)
                top_score = m.get("results", {}).get("top_score")
                val_stat = m.get("status", val_stat)
            except Exception:
                pass

        runs.append({
            "run_id": d.name,
            "dir": d,
            "target_id": target_id,
            "ligand_name": lig_name,
            "top_score": top_score,
            "status": val_stat,
        })
    return runs


def scan_completed_md_runs(root: Path) -> List[Path]:
    """Return MD run directories whose manifest records a completed run."""
    if not root.is_dir():
        return []
    found: List[Path] = []
    for directory in sorted(root.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True):
        if not directory.is_dir():
            continue
        manifest_file = directory / "md_manifest.json"
        if not manifest_file.is_file():
            continue
        try:
            manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if manifest.get("status") == "completed":
            found.append(directory)
    return found


def load_saved_md_result(run_dir: Path) -> Optional[Dict[str, Any]]:
    """Rebuild the session payload for a completed run stored on disk."""
    manifest_file = run_dir / "md_manifest.json"
    try:
        manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if manifest.get("status") != "completed":
        return None
    topology = manifest.get("topology") or {}
    pose_meta = {
        "ligand_name": topology.get("ligand_resname") or "Ligand",
        "target_pdb": manifest.get("run_id") or run_dir.name,
    }
    return {
        "run_id": manifest.get("run_id") or run_dir.name,
        "run_dir": run_dir,
        "manifest": manifest,
        "pose_meta": pose_meta,
        "ligand_resname": topology.get("ligand_resname") or "LIG",
    }


def _derive_ligand_resname(ligand_name: str, sdf_path: Optional[str]) -> str:
    if sdf_path:
        stem = Path(sdf_path).stem
        token = re.match(r"[A-Za-z0-9]{1,5}", stem)
        if token:
            return sanitize_resname(token.group(0))
    token = re.match(r"[A-Za-z0-9]{1,5}", ligand_name or "")
    return sanitize_resname(token.group(0) if token else "LIG")


def _write_failure_manifest(run_dir: Path, error: str) -> None:
    manifest = {
        "run_id": run_dir.name,
        "status": "failed",
        "error": error,
        "protocol_reference": "wiki/gromacs-protein-ligand-protocol.md",
    }
    try:
        (run_dir / "md_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
        )
    except OSError:
        pass


def _build_mdp_bundle(target_dir: Path, cfg: SimulationConfig) -> None:
    mdp = {
        "ions.mdp": generate_mdp_ions(),
        "em.mdp": generate_mdp_minimization(cfg.force_field, cfg.em_tol, cfg.em_max_steps),
        "nve.mdp": generate_mdp_nve(
            cfg.nve_duration_ps, cfg.temperature_k, cfg.timestep_ps,
            cfg.force_field, cfg.random_seed,
        ),
        "nvt.mdp": generate_mdp_nvt(
            cfg.temperature_k, cfg.nvt_duration_ps, cfg.timestep_ps, cfg.force_field, cfg.random_seed
        ),
        "npt.mdp": generate_mdp_npt(
            cfg.temperature_k, cfg.pressure_bar, cfg.npt_duration_ps, cfg.timestep_ps, cfg.force_field
        ),
        "md.mdp": generate_mdp_production(
            cfg.prod_duration_ns, cfg.temperature_k, cfg.pressure_bar, cfg.timestep_ps, cfg.force_field
        ),
    }
    for replicate_index in range(2, max(1, min(int(cfg.n_replicates), 5)) + 1):
        mdp[f"nvt_rep{replicate_index}.mdp"] = generate_mdp_nvt(
            temperature_k=cfg.temperature_k,
            duration_ps=cfg.nvt_duration_ps,
            timestep_ps=cfg.timestep_ps,
            force_field=cfg.force_field,
            random_seed=cfg.random_seed + replicate_index * 1000,
        )
    for name, content in mdp.items():
        (target_dir / name).write_text(content, encoding="utf-8", newline="\n")


def render_md_tab():
    """Main entrypoint for Tab 3: Molecular Dynamics Simulation & Analysis."""
    st.header("🔬 MD Simulation — GROMACS Pipeline")

    gmx_info = get_gmx_version()
    backend = detect_gromacs()
    col_stat1, col_stat2, col_stat3 = st.columns(3)
    with col_stat1:
        if gmx_info["available"]:
            st.success(f"🟢 GROMACS Available ({gmx_info['version']})")
        else:
            st.error("🔴 GROMACS unavailable. No simulation will be executed.")
    with col_stat2:
        st.metric(
            "CUDA Acceleration",
            "Enabled" if gmx_info.get("gpu_acceleration") else "CPU / SIMD",
        )
    with col_stat3:
        st.metric("Unit Cell", "Rhombic Dodecahedron")

    if backend is None:
        st.info(
            "💡 Install GROMACS natively or inside "
            "WSL2 Ubuntu and verify with `gmx --version`. "
            "You may also configure the `GMX_BIN` environment variable."
        )
        return

    tools = check_backend_tools(backend)
    if not tools.get("acpype"):
        st.error(
            "🔴 Small-molecule parameterization unavailable: ACPYPE "
            "is not installed in the backend environment. Install it inside WSL2 "
            "Ubuntu (e.g., `pip install --user acpype`); ACPYPE bundles its own AmberTools binaries."
        )
        return

    st.markdown("---")
    st.subheader("1️⃣ Docking Run Selection & Optimal Binding Pose Extraction")

    current_docking = st.session_state.get("current_docking_result")
    saved_runs = scan_available_docking_runs()

    run_options: List[str] = []
    run_dir_map: Dict[str, Path] = {}
    run_labels: Dict[str, str] = {}
    if current_docking and "run_dir" in current_docking:
        cid = current_docking.get("run_id", "active_session")
        opt_key = f"active::{cid}" if any(r["run_id"] == cid for r in saved_runs) else str(cid)
        lig = current_docking.get("ligand_name") or current_docking.get("ligand", {}).get("ligand_name") or "Ligand"
        pdb = current_docking.get("target_id") or current_docking.get("target", {}).get("pdb_id") or "Target"
        score = current_docking.get("top_score")
        if score is None:
            score = current_docking.get("best_vina_score")
        if score is None and current_docking.get("poses"):
            first_pose = current_docking["poses"][0]
            score = getattr(first_pose, "score", None) or (first_pose.get("score") if isinstance(first_pose, dict) else None)
        score_part = f" ({score:.2f} kcal/mol)" if score is not None else ""
        label = f"Active Session — {lig} × {pdb}{score_part}"
        run_options.append(opt_key)
        run_dir_map[opt_key] = Path(current_docking["run_dir"])
        run_labels[opt_key] = label

    for r in saved_runs:
        opt_key = str(r["run_id"])
        if opt_key in run_dir_map:
            continue
        score_part = f" [{r['top_score']:.2f} kcal/mol]" if r.get("top_score") is not None else ""
        match = re.search(r"(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})", opt_key)
        if match:
            ts = f"{match.group(1)}-{match.group(2)}-{match.group(3)} {match.group(4)}:{match.group(5)}"
        else:
            try:
                ts = datetime.fromtimestamp(r["dir"].stat().st_mtime).strftime("%Y-%m-%d %H:%M")
            except Exception:
                ts = ""
        ts_part = f" · {ts}" if ts else ""

        lower_id = opt_key.lower()
        if "matched_mutant" in lower_id:
            label = f"Mutant ({r['target_id']}) — {r['ligand_name']}{score_part}{ts_part}"
        elif "matched_healthy" in lower_id:
            label = f"Healthy ({r['target_id']}) — {r['ligand_name']}{score_part}{ts_part}"
        else:
            label = f"{r['ligand_name']} × {r['target_id']}{score_part}{ts_part}"

        run_options.append(opt_key)
        run_dir_map[opt_key] = r["dir"]
        run_labels[opt_key] = label

    if not run_options:
        st.warning("⚠️ No completed docking runs were found in `docking_runs/` or in the active session.")
        st.info("💡 Run molecular docking first in Tab 2 ('💊 Tab 2: Molecular Docking Simulation (AutoDock Vina)'), then return here. No synthetic data will be generated.")
        return

    selected_run_id = st.selectbox(
        "Select Source Docking Run:",
        run_options,
        format_func=lambda opt: run_labels.get(opt, str(opt)),
    )
    selected_run_dir = run_dir_map[selected_run_id]

    try:
        rec_path, best_pose, pose_meta = extract_best_pose_from_docking_run(selected_run_dir)
    except Exception as ex:  # noqa: BLE001
        st.error(f"Error extracting binding pose from docking run: {ex}")
        return

    col_meta1, col_meta2, col_meta3, col_meta4 = st.columns(4)
    with col_meta1:
        st.metric("Target PDB", pose_meta.get("target_pdb", "Unknown"))
    with col_meta2:
        st.metric("Ligand Name", pose_meta.get("ligand_name", "Ligand"))
    with col_meta3:
        st.metric("Best Docking Score", f"{pose_meta['best_vina_score']:.2f} kcal/mol")
    with col_meta4:
        st.metric("Heavy Atom Count", f"{pose_meta['heavy_atoms']} (total {pose_meta['total_atoms']})")

    if not pose_meta.get("ligand_sdf_path"):
        st.warning(
            "⚠️ No ligand SDF was found in the docking run directory. "
            "The pipeline will attempt to infer bond orders from the pose, but re-running docking produces a more accurate topology."
        )

    with st.expander("👁️ Preview Extracted Binding Pose (3D Complex)", expanded=False):
        rec_text = Path(rec_path).read_text(encoding="utf-8", errors="replace")
        combined_preview_pdb = rec_text.replace("END\n", "") + "\n" + best_pose.pdb_block
        render_mol_3d(combined_preview_pdb, height=360, ligand_resname="LIG")

    st.markdown("---")
    st.subheader("2️⃣ Complex Topology & Force Field Configuration")

    amber_fields = [ff for ff in tools.get("force_fields", []) if ff.startswith("amber")]
    if not amber_fields:
        st.error(
            "No AMBER force field is installed "
            "for this backend. Install `amber99sb-ildn.ff` or use a full GROMACS installation."
        )
        return
    default_ff = "amber99sb-ildn" if "amber99sb-ildn" in amber_fields else amber_fields[0]

    col_top1, col_top2 = st.columns(2)
    with col_top1:
        ff_choice = st.selectbox(
            "Protein Force Field:",
            amber_fields,
            index=amber_fields.index(default_ff),
            help="Supported route: Amber99SB-ILDN for the protein + GAFF2 / AM1-BCC for the ligand.",
        )
    with col_top2:
        water_options = {
            "TIP3P": "tip3p",
            "SPC/E": "spc",
            "TIP4P-Ew": "tip4pew",
        }
        water_label = st.selectbox(
            "Explicit Water Model:",
            list(water_options.keys()),
        )
    water_model = water_options[water_label]

    st.markdown("---")
    st.subheader("3️⃣ Simulation Box & Simulation Parameters")

    col_box1, col_box2 = st.columns(2)
    with col_box1:
        box_label = st.selectbox(
            "Unit Cell Geometry:",
            ["dodecahedron (Rhombic Dodecahedron)", "cubic (Cubic Box)"],
        )
    with col_box2:
        box_dist = st.slider("Solvent Padding Distance (nm):", 0.8, 1.5, 1.0, 0.1)

    col_sim1, col_sim2, col_sim3 = st.columns(3)
    with col_sim1:
        temp_k = st.number_input("Temperature (K):", value=300.0, step=5.0)
    with col_sim2:
        press_bar = st.number_input("Pressure (bar):", value=1.0, step=0.1)
    with col_sim3:
        prod_ns = st.slider("Production Duration (ns):", 0.01, 10.0, 0.5, 0.01)

    if prod_ns < 1.0:
        st.warning(
            "⚠️ Production duration is below 1.0 ns. Trajectories this short are "
            "generally insufficient for equilibrated stability or convergence claims; "
            "treat RMSD, RMSF, and H-bond results as exploratory only."
        )

    col_nve1, col_nve2 = st.columns([1.0, 1.0])
    with col_nve1:
        nve_enabled = st.checkbox(
            "Run NVE drift check (LiveCoMS)",
            value=True,
            help="Short uncoupled run without thermostat or barostat to verify absence of systematic energy drift.",
        )
    with col_nve2:
        nve_ps = st.number_input(
            "NVE Check Duration (ps):",
            min_value=1.0,
            max_value=100.0,
            value=10.0,
            step=1.0,
            disabled=not nve_enabled,
        )
    n_replicates = st.slider(
        "Independent Production Replicates:",
        min_value=1,
        max_value=5,
        value=3,
        help=(
            "Replicate 1 is the primary trajectory. Additional replicates start from the same "
            "EM structure with independent Maxwell–Boltzmann velocities and repeat "
            "NVT, NPT, and production."
        ),
    )
    allow_cofactor_removal = st.checkbox(
        "Allow removal of metallic cofactors (HEM/FE/MG/ZN/MN/CU/CA)",
        value=False,
        help=(
            "By default, MD is refused if essential cofactors are present because "
            "removing them alters system chemistry. Enable this only after confirming the cofactor is non-essential."
        ),
    )
    retain_monatomic_ions = st.checkbox(
        "Retain monatomic metal ions (ZN/MG/CA/K)",
        value=False,
        help=(
            "When enabled, supported monatomic ions are retained in the simulated system instead of being stripped. "
            "Polyatomic cofactors such as HEM remain rejected."
        ),
    )
    retain_crystal_waters = st.checkbox(
        "Retain crystallographic water molecules (HOH/WAT)",
        value=False,
        help=(
            "Crystallographic waters are retained as solvent residues. "
            "Missing hydrogens are placed geometrically using TIP3P geometry with approximate orientation."
        ),
    )

    st.markdown(
        "**Thermostat:** `v-rescale` | "
        "**Coupling Groups (`tc-grps`):** `Protein_Ligand Water_and_ions`"
    )
    st.markdown(
        "**NPT Equilibration:** `berendsen` — 100 ps | "
        "**Production Barostat:** `Parrinello-Rahman` | "
        "**Electrostatics:** PME | "
        "**Bond Constraints:** LINCS"
    )

    ligand_resname = _derive_ligand_resname(
        pose_meta.get("ligand_name", "Ligand"), pose_meta.get("ligand_sdf_path")
    )
    cfg = SimulationConfig(
        run_id="",
        temperature_k=float(temp_k),
        pressure_bar=float(press_bar),
        prod_duration_ns=float(prod_ns),
        nve_duration_ps=float(nve_ps) if nve_enabled else 0.0,
        n_replicates=int(n_replicates),
        box_type="dodecahedron" if "dodecahedron" in box_label else "cubic",
        box_distance_nm=float(box_dist),
        force_field=ff_choice,
        water_model=water_model,
    )

    st.markdown("---")
    col_btn1, col_btn2 = st.columns(2)
    run_clicked = col_btn1.button(
        "🚀 Run MD Pipeline",
        type="primary",
        width="stretch",
    )
    export_clicked = col_btn2.button(
        "📥 Export Complete GROMACS Protocol Bundle",
        width="stretch",
    )

    def _build_topology(target_dir: Path):
        sdf = pose_meta.get("ligand_sdf_path")
        return prepare_complex_topology(
            receptor_pdb=rec_path,
            ligand_pose=best_pose,
            output_dir=target_dir,
            force_field=ff_choice,
            water_model=water_model,
            ligand_sdf=Path(sdf) if sdf else None,
            ligand_resname=ligand_resname,
            backend=backend,
            allow_cofactor_removal=allow_cofactor_removal,
            retain_monatomic_ions=retain_monatomic_ions,
            retain_crystal_waters=retain_crystal_waters,
        )

    if export_clicked:
        with tempfile.TemporaryDirectory() as tmp_export:
            tmp_path = Path(tmp_export)
            try:
                with st.spinner("Building real topology via pdb2gmx and ACPYPE..."):
                    topology = _build_topology(tmp_path)
                    _build_mdp_bundle(tmp_path, cfg)
                    generate_standalone_scripts(
                        output_dir=tmp_path,
                        complex_gro=topology.complex_gro.name,
                        topol_top=topology.topol_top.name,
                        index_ndx="index.ndx",
                        box_type=cfg.box_type,
                        box_distance=cfg.box_distance_nm,
                        production_prefix="prod_md",
                        ligand_resname=ligand_resname,
                        water_model=water_model,
                        n_replicates=cfg.n_replicates,
                        salt_concentration_molar=cfg.salt_concentration_molar,
                    )
                    shutil.copy2(Path(__file__).resolve().parent / "indexing.py", tmp_path / "make_index.py")
                    shutil.copy2(Path(__file__).resolve().parent / "finalize_run.py", tmp_path / "finalize_md_run.py")
                    (tmp_path / "README_BUNDLE.txt").write_text(
                        "Bio_V2 GROMACS reproducibility bundle\n"
                        "====================================\n"
                        "index.ndx is NOT shipped as a static file: it is generated from the actual\n"
                        "solv_ions.gro by make_index.py in step 4 of run_md.sh / run_md.bat.\n"
                        "This prevents atom-count mismatches between the index and the real system.\n\n"
                        "Run ./run_md.sh (Linux/WSL) or run_md.bat (native Windows GROMACS).\n",
                        encoding="utf-8",
                        newline="\n",
                    )

                zip_target = tmp_path / "gromacs_md_bundle.zip"
                with zipfile.ZipFile(zip_target, "w", zipfile.ZIP_DEFLATED) as zipf:
                    for p in tmp_path.rglob("*"):
                        if p.is_file() and p.suffix != ".zip":
                            zipf.write(p, arcname=str(p.relative_to(tmp_path)))
                zip_bytes = zip_target.read_bytes()

                st.download_button(
                    label="⬇️ Download GROMACS Bundle (.ZIP)",
                    data=zip_bytes,
                    file_name=f"md_bundle_{pose_meta.get('target_pdb', 'target')}_{ligand_resname}.zip",
                    mime="application/zip",
                    width="stretch",
                )
                st.success("✅ Real protocol bundle ready for download.")
                if topology.warnings:
                    for w in topology.warnings:
                        st.warning(w)
            except (TopologyBuildError, Exception) as exc:  # noqa: BLE001
                st.error(f"Failed to build GROMACS bundle: {exc}")

    if run_clicked:
        run_id, md_run_dir = create_md_run_directory(prefix="md")
        progress_bar = st.progress(0.0)
        status_text = st.empty()

        def on_progress(stage: str, frac: float):
            progress_bar.progress(min(max(frac, 0.0), 1.0))
            status_text.markdown(f"**Current Stage:** {stage}")

        cfg.run_id = run_id
        try:
            with st.spinner("Building real topology (pdb2gmx + ACPYPE)..."):
                topology = _build_topology(md_run_dir)
            engine = BiophysicalMDSimulator(config=cfg, output_dir=md_run_dir, backend=backend)
            manifest = engine.run_pipeline(topology=topology, progress_callback=on_progress)
            st.session_state["current_md_result"] = {
                "run_id": run_id,
                "run_dir": md_run_dir,
                "manifest": manifest,
                "complex_pdb": topology.complex_pdb,
                "topology": topology,
                "best_pose": best_pose,
                "pose_meta": pose_meta,
                "config": cfg,
                "ligand_resname": ligand_resname,
            }
            st.success(f"✅ Real GROMACS simulation completed successfully! ({run_id})")
        except (TopologyBuildError, SimulationError, Exception) as exc:  # noqa: BLE001
            _write_failure_manifest(md_run_dir, str(exc))
            st.session_state.pop("current_md_result", None)
            st.error("❌ Simulation failed. No synthetic results were generated.")
            st.exception(exc)

    # ------------------------------------------------------------------
    # Results (only real GROMACS artifacts)
    # ------------------------------------------------------------------
    saved_runs = scan_completed_md_runs(MD_RUNS_DIR)
    if saved_runs:
        choices = ["Current session"] + [directory.name for directory in saved_runs]
        selected = st.selectbox(
            "Load a saved MD run (completed only):",
            choices,
            key="md_saved_run_choice",
        )
        if selected != "Current session":
            md_res = load_saved_md_result(MD_RUNS_DIR / selected)
        else:
            md_res = st.session_state.get("current_md_result")
    else:
        md_res = st.session_state.get("current_md_result")
    if not md_res:
        return

    manifest = md_res.get("manifest", {})
    if manifest.get("status") != "completed":
        st.error("This run did not complete successfully; no results are available to display.")
        return

    st.markdown("---")
    st.subheader(
        f"📊 Real Trajectory Results: {md_res['pose_meta'].get('ligand_name')} ⚯ "
        f"{md_res['pose_meta'].get('target_pdb')}"
    )

    run_dir = Path(md_res["run_dir"])
    lig_name = md_res["pose_meta"].get("ligand_name", "Ligand")

    # RMSD / RMSF from real gmx rms / gmx rmsf XVG (nm -> Angstrom)
    rmsd_records: List[RMSDRecord] = []
    rmsd_xvg = run_dir / "rmsd_ligand.xvg"
    for t, val in parse_xvg_file(rmsd_xvg):
        rmsd_records.append(RMSDRecord(time_ns=t, rmsd_angstrom=val * 10.0))

    rmsf_records: List[RMSFRecord] = []
    rmsf_xvg = run_dir / "rmsf_protein.xvg"
    for res_idx, val in parse_xvg_file(rmsf_xvg):
        rmsf_records.append(
            RMSFRecord(resnum=int(res_idx), resname="RES", rmsf_angstrom=val * 10.0)
        )

    hbond_points = parse_hbond_num_xvg(run_dir / "hbond_num.xvg")
    mmpbsa_res = load_mmpbsa_results(run_dir)

    mean_rmsd = None
    mean_hbonds = None
    burn_in_fraction = float(
        manifest.get("parameters", {}).get(
            "analysis_burn_in_fraction",
            getattr(md_res.get("config"), "analysis_burn_in_fraction", 0.2),
        )
    )
    rmsd_summary = summarize_timeseries(
        [(r.time_ns, r.rmsd_angstrom) for r in rmsd_records],
        burn_in_fraction=burn_in_fraction,
    )
    hbond_summary = summarize_timeseries(hbond_points, burn_in_fraction=burn_in_fraction)
    if rmsd_summary["n"]:
        mean_rmsd = rmsd_summary["mean"]
    if hbond_summary["n"]:
        mean_hbonds = hbond_summary["mean"]

    col_kpi1, col_kpi2, col_kpi3, col_kpi4 = st.columns(4)
    with col_kpi1:
        st.metric(
            "Mean Ligand RMSD",
            f"{mean_rmsd:.2f} Å" if mean_rmsd is not None else "Not available",
            ("Stable (≤ 2.0 Å)" if mean_rmsd <= 2.0 else "Notable Fluctuation")
            if mean_rmsd is not None else None,
        )
    with col_kpi2:
        st.metric(
            "Mean H-bond Count",
            f"{mean_hbonds:.1f}" if mean_hbonds is not None else "Not available",
        )
    with col_kpi3:
        production_ns = manifest.get("parameters", {}).get("production_duration_ns")
        st.metric(
            "Production Duration",
            f"{production_ns:.2f} ns" if production_ns is not None else "Not recorded",
        )
    with col_kpi4:
        engine_ver = manifest.get("engine", {}).get("version", "").strip()
        engine_label = f"GROMACS {engine_ver}" if engine_ver and engine_ver != "unknown" else "GROMACS 2026.3"
        st.metric(
            engine_label,
            "Completed" if manifest.get("status") == "completed" else "Ready",
        )
    uncertainty = []
    for label, summary, unit in (("RMSD", rmsd_summary, "Å"), ("H-bonds", hbond_summary, "")):
        if summary["n"]:
            uncertainty.append(
                f"{label} = {summary['mean']:.2f} ± {summary['sem']:.2f} {unit} "
                f"(block SEM {summary['block_sem']:.2f}, n={int(summary['n'])})"
            )
    if uncertainty:
        st.caption(
            f"Within-trajectory uncertainty (burn-in {burn_in_fraction * 100:.0f}%): "
            + "; ".join(uncertainty)
            + ". Not inter-replicate variance."
        )
    replicate_stats = manifest.get("replicate_statistics") or {}
    if int(replicate_stats.get("n_replicates", 1)) > 1:
        rmsd_rep = replicate_stats.get("rmsd_angstrom", {})
        hb_rep = replicate_stats.get("hbond_count", {})
        comparisons = []
        for label, summary, unit in (("RMSD", rmsd_rep, "Å"), ("H-bonds", hb_rep, "")):
            if summary.get("n", 0) and summary.get("mean") is not None:
                comparisons.append(
                    f"{label} = {summary['mean']:.2f} ± {summary['sd']:.2f} {unit} "
                    f"(n={int(summary['n'])})"
                )
        if comparisons:
            st.caption("Across independent replicates: " + " | ".join(comparisons) + ".")

    tab_rmsd, tab_rmsf, tab_hbond, tab_mmpbsa, tab_traj3d = st.tabs([
        "📈 Pose Stability — RMSD",
        "🧬 Protein Flexibility — RMSF",
        "🔗 Hydrogen Bonds",
        "⚡ Binding Energy — MM-PB(GB)SA",
        "🎥 Trajectory Viewer (3D)",
    ])

    with tab_rmsd:
        st.plotly_chart(create_rmsd_plotly(rmsd_records, lig_name), width="stretch")
        st.caption("Reference: ≤ 2.0 Å stable; > 3.5 Å drift.")

    with tab_rmsf:
        st.plotly_chart(create_rmsf_plotly(rmsf_records, md_res['pose_meta'].get('target_pdb')), width="stretch")

    with tab_hbond:
        st.plotly_chart(create_hbond_count_plotly(hbond_points, lig_name), width="stretch")
        st.caption("Source: gmx hbond -num (PBC-corrected; cutoff: 0.35 nm, 30°).")

    with tab_mmpbsa:
        if mmpbsa_res is None:
            st.warning(MMPBSA_NOT_AVAILABLE_MESSAGE)
            st.caption(
                "Run gmx_MMPBSA manually on the PBC-wrapped production trajectory to compute binding free energy."
            )
        else:
            st.plotly_chart(create_mmpbsa_breakdown_plotly(mmpbsa_res), width="stretch")
            st.markdown(
                f"**Solvation Model:** {mmpbsa_res.solvation_model} | "
                f"**ΔG_bind:** {mmpbsa_res.delta_g_bind:.2f} ± {mmpbsa_res.delta_g_bind_std:.2f} kcal/mol"
            )
            eps_text = (
                f"{mmpbsa_res.internal_dielectric:g}"
                if mmpbsa_res.internal_dielectric_specified
                else "Not specified"
            )
            st.caption(
                f"Frames evaluated: {mmpbsa_res.frames_evaluated or 'Not specified'} | "
                f"eps_int (internal dielectric): {eps_text} | "
                f"GB model: {mmpbsa_res.igb if mmpbsa_res.igb is not None else 'Not specified'} | "
                f"-TΔS (entropy term, kcal/mol): {mmpbsa_res.entropy_term_kcal_mol if mmpbsa_res.entropy_term_kcal_mol is not None else 'Not specified'}"
            )

    with tab_traj3d:
        st.markdown("**Final Production Frame — PBC-wrapped & fitted:**")
        last_frame = load_last_frame_pdb(run_dir / "trajectory_pbc_fit.pdb")
        if last_frame:
            c_surf1, c_surf2, c_surf3 = st.columns([1, 1, 1])
            with c_surf1:
                show_surface = st.checkbox("Show Receptor Surface", value=True, key="md_show_surface")
            if show_surface:
                with c_surf2:
                    surface_type = st.selectbox(
                        "Surface Model",
                        ["MS (Molecular Surface)", "SAS (Solvent-Accessible Surface)", "VDW (Van der Waals Surface)"],
                        key="md_surface_type",
                    )
                with c_surf3:
                    surface_opacity = st.slider(
                        "Surface Opacity",
                        min_value=0.1,
                        max_value=1.0,
                        value=0.5,
                        step=0.05,
                        key="md_surface_opacity",
                    )
            else:
                surface_type = "MS"
                surface_opacity = 0.5

            render_mol_3d(
                last_frame,
                height=450,
                ligand_resname=md_res.get("ligand_resname", "LIG"),
                show_surface=show_surface,
                surface_opacity=surface_opacity,
                surface_type=surface_type.split()[0] if isinstance(surface_type, str) else "MS",
            )
        else:
            st.info("Trajectory PDB file not found.")

    st.caption(f"🛡️ {DISCLAIMER_MD}")
