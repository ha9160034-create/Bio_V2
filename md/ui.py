"""
Molecular Dynamics Streamlit UI & Interactive Biophysical Dashboard
===================================================================

Scientific Traceability & Wiki Documentation:
---------------------------------------------
- `wiki/gromacs-protein-ligand-protocol.md`: Justin Lemkul's GROMACS Protocol.
- `wiki/open-force-field-and-ligand-parameterization.md`: Small-Molecule Force Fields (OpenFF Sage, CGenFF, ACPYPE).
- `wiki/simulation-and-benchmark-best-practices.md`: LiveCoMS Foundations & Benchmark Best Practices.
- `wiki/mm-pbsa-and-mm-gbsa.md`: MM-PBSA / MM-GBSA Binding Free Energy Rescoring.

Biophysical Implementation Guardrails:
--------------------------------------
1. Unit Cell: Default to rhombic dodecahedron (-bt dodecahedron -d 1.0) to save ~23% solvent volume.
2. Restraints: Harmonic position restraints (k = 1000 kJ/(mol*nm^2)) on ligand heavy atoms during equilibration.
3. Thermostats: Merged tc-grps = Protein_Ligand Water_and_ions (never couple small ligand independently).
4. Ensemble Validity: Reject Berendsen barostat in production; enforce Parrinello-Rahman and v-rescale.
5. Trajectory Processing: Always apply PBC re-wrapping (trjconv -pbc mol -ur compact) before analysis.

Honesty policy
--------------
This tab only reports artifacts produced by a real GROMACS run. When GROMACS,
ACPYPE or AmberTools are missing, the UI explains what is missing and refuses
to run instead of fabricating a trajectory or an MM-PBSA number.
"""

from __future__ import annotations

import re
import json
import shutil
import zipfile
import tempfile
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
    "Biophysical Scientific Notice: MD (Molecular Dynamics) results depend on classical "
    "force fields and implicit approximations for ions and solvent, whereas MM-PB(GB)SA — "
    "MM-PBSA (Molecular Mechanics / Poisson–Boltzmann Surface Area) / "
    "MM-GBSA (Molecular Mechanics / Generalized Born Surface Area) — values are displayed only "
    "when gmx_MMPBSA has actually been executed. Results are intended for stability evaluation "
    "and relative compound ranking, not as a substitute for experimental assays "
    "(ITC — Isothermal Titration Calorimetry / SPR — Surface Plasmon Resonance)."
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
                lig_name = m.get("ligand", {}).get("name", lig_name)
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
    st.header("🔬 MD (Molecular Dynamics) — GROMACS (GROningen MAchine for Chemical Simulations) Pipeline")
    st.caption("All-atom protein–ligand simulation protocol following Lemkul and LiveCoMS standards documented in wiki/")

    with st.expander("📖 Scientific References & Protocol Standards (Wiki & Biophysical Guardrails)", expanded=False):
        st.markdown("""
        - **`wiki/gromacs-protein-ligand-protocol.md`**: Lemkul protocol for protein–ligand complex preparation, equilibration, and production MD (Molecular Dynamics).
        - **`wiki/open-force-field-and-ligand-parameterization.md`**: Small-molecule parameterization via ACPYPE (AnteChamber PYthon Parser interfacE), GAFF2 (General AMBER Force Field v2), and AM1-BCC (Austin Model 1 with Bond Charge Corrections).
        - **`wiki/simulation-and-benchmark-best-practices.md`**: Prohibition of the Berendsen barostat in production and merged thermostat coupling groups.
        - **`wiki/mm-pbsa-and-mm-gbsa.md`**: MM-PB(GB)SA — MM-PBSA (Molecular Mechanics / Poisson–Boltzmann Surface Area) / MM-GBSA (Molecular Mechanics / Generalized Born Surface Area) — computed exclusively via gmx_MMPBSA on the real production trajectory.
        """)

    gmx_info = get_gmx_version()
    backend = detect_gromacs()
    col_stat1, col_stat2, col_stat3 = st.columns(3)
    with col_stat1:
        if gmx_info["available"]:
            where = "WSL2 (Windows Subsystem for Linux v2)" if gmx_info.get("backend") == "wsl" else "Native"
            st.success(f"🟢 GROMACS (GROningen MAchine for Chemical Simulations) Available ({gmx_info['version']}) — {where}")
        else:
            st.error("🔴 GROMACS (GROningen MAchine for Chemical Simulations) unavailable. No simulation will be executed.")
    with col_stat2:
        st.metric(
            "Hardware Acceleration — CUDA (Compute Unified Device Architecture GPU Acceleration)",
            "Enabled" if gmx_info.get("gpu_acceleration") else "CPU / SIMD (Single Instruction, Multiple Data CPU Vectorization)",
        )
    with col_stat3:
        st.metric("Default Unit Cell Geometry", "Rhombic Dodecahedron (-23% solvent)")

    if backend is None:
        st.info(
            "💡 Install GROMACS (GROningen MAchine for Chemical Simulations) natively or inside "
            "WSL2 (Windows Subsystem for Linux v2) Ubuntu and verify with `gmx --version`. "
            "You may also configure the `GMX_BIN` environment variable."
        )
        return

    tools = check_backend_tools(backend)
    if not tools.get("acpype"):
        st.error(
            "🔴 Small-molecule parameterization unavailable: ACPYPE (AnteChamber PYthon Parser interfacE) "
            "is not installed in the backend environment. Install it inside WSL2 (Windows Subsystem for Linux v2) "
            "Ubuntu (e.g., `pip install --user acpype`); ACPYPE bundles its own AmberTools binaries."
        )
        return

    with st.expander("🧰 Simulation Backend Preflight Check", expanded=False):
        st.markdown(
            f"- **GROMACS (GROningen MAchine for Chemical Simulations):** `{gmx_info['path']}` ({gmx_info.get('backend')})\n"
            f"- **ACPYPE (AnteChamber PYthon Parser interfacE):** {'✅' if tools.get('acpype') else '❌'} | "
            f"**Antechamber:** {'✅ (PATH or bundled ACPYPE)' if tools.get('antechamber') else '❌'} | "
            f"**OpenBabel:** {'✅' if tools.get('obabel') else '❌'}\n"
            f"- **Available Force Fields:** {', '.join(tools.get('force_fields', []))}"
        )

    st.markdown("---")
    st.subheader("1️⃣ Docking Run Selection & Optimal Binding Pose Extraction")

    current_docking = st.session_state.get("current_docking_result")
    saved_runs = scan_available_docking_runs()

    run_options: List[str] = []
    run_dir_map: Dict[str, Path] = {}
    if current_docking and "run_dir" in current_docking:
        cid = current_docking.get("run_id", "active_session")
        label = (
            f"Active Session ({cid}) — {current_docking.get('ligand_name', 'Ligand')} "
            f"⚯ {current_docking.get('target_id', 'Target')}"
        )
        run_options.append(label)
        run_dir_map[label] = Path(current_docking["run_dir"])

    for r in saved_runs:
        score_str = f"[{r['top_score']:.2f} kcal/mol]" if r["top_score"] is not None else ""
        label = f"{r['run_id']} — {r['ligand_name']} ⚯ {r['target_id']} {score_str}"
        if label not in run_dir_map:
            run_options.append(label)
            run_dir_map[label] = r["dir"]

    if not run_options:
        st.warning("⚠️ No completed docking runs were found in `docking_runs/` or in the active session.")
        st.info("💡 Run molecular docking first in Tab 2 ('💊 Tab 2: Molecular Docking Simulation (AutoDock Vina)'), then return here. No synthetic data will be generated.")
        return

    selected_run_label = st.selectbox("Select Source Docking Run:", run_options)
    selected_run_dir = run_dir_map[selected_run_label]

    try:
        rec_path, best_pose, pose_meta = extract_best_pose_from_docking_run(selected_run_dir)
    except Exception as ex:  # noqa: BLE001
        st.error(f"Error extracting binding pose from docking run: {ex}")
        return

    col_meta1, col_meta2, col_meta3, col_meta4 = st.columns(4)
    with col_meta1:
        st.metric("Target — PDB (Protein Data Bank)", pose_meta.get("target_pdb", "Unknown"))
    with col_meta2:
        st.metric("Ligand Name", pose_meta.get("ligand_name", "Ligand"))
    with col_meta3:
        st.metric("Best Docking Score — kcal/mol (kilocalories per mole)", f"{pose_meta['best_vina_score']:.2f} kcal/mol")
    with col_meta4:
        st.metric("Heavy Atom Count", f"{pose_meta['heavy_atoms']} (total {pose_meta['total_atoms']})")

    if not pose_meta.get("ligand_sdf_path"):
        st.warning(
            "⚠️ No ligand SDF (Structure-Data File) was found in the docking run directory. "
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
            "No AMBER (Assisted Model Building with Energy Refinement Force Field) force field is installed "
            "for this backend. Install `amber99sb-ildn.ff` or use a full GROMACS installation."
        )
        return
    default_ff = "amber99sb-ildn" if "amber99sb-ildn" in amber_fields else amber_fields[0]

    col_top1, col_top2 = st.columns(2)
    with col_top1:
        ff_choice = st.selectbox(
            "Protein Force Field — AMBER (Assisted Model Building with Energy Refinement Force Field):",
            amber_fields,
            index=amber_fields.index(default_ff),
            help=(
                "Supported route: Amber99SB-ILDN for the protein + "
                "GAFF2 (General AMBER Force Field v2) / AM1-BCC (Austin Model 1 with Bond Charge Corrections) for the ligand."
            ),
        )
    with col_top2:
        water_options = {
            "TIP3P (Transferable Intermolecular Potential with 3 Points Water Model)": "tip3p",
            "SPC/E (Extended Simple Point Charge Water Model)": "spc",
            "TIP4P-Ew (Transferable Intermolecular Potential with 4 Points Ewald Water Model)": "tip4pew",
        }
        water_label = st.selectbox(
            "Explicit Water Model:",
            list(water_options.keys()),
        )
    water_model = water_options[water_label]

    st.info(
        "🔒 **Only Real Parameterization Route Available:** `gmx pdb2gmx` for the protein + "
        "ACPYPE (AnteChamber PYthon Parser interfacE) using GAFF2 (General AMBER Force Field v2) + "
        "AM1-BCC (Austin Model 1 with Bond Charge Corrections) for the ligand. "
        "OpenFF (Open Force Field Sage) and CHARMM36 (Chemistry at HARvard Macromolecular Mechanics Force Field v36) / "
        "CGenFF (CHARMM General Force Field) options are disabled because their tools/files are not installed, "
        "and placeholder parameters will never be generated."
    )
    st.info(
        "🔒 **Heavy-Atom Position Restraints:** Harmonic restraints with force constant "
        "1000 kJ/(mol·nm²) (kilojoules per mole per square nanometer) are applied to ligand heavy atoms "
        "during NVT (Canonical Ensemble: Constant Particles, Volume, and Temperature) and "
        "NPT (Isothermal–Isobaric Ensemble: Constant Particles, Pressure, and Temperature) equilibration "
        "to prevent pocket drift."
    )

    st.markdown("---")
    st.subheader("3️⃣ Simulation Box & MD (Molecular Dynamics) Parameters")

    col_box1, col_box2 = st.columns(2)
    with col_box1:
        box_label = st.selectbox(
            "Unit Cell Geometry:",
            ["dodecahedron (Rhombic Dodecahedron — Recommended)", "cubic (Cubic Box)"],
        )
    with col_box2:
        box_dist = st.slider("Solvent Padding Distance, nm (nanometers):", 0.8, 1.5, 1.0, 0.1)

    if "dodecahedron" in box_label:
        _, _, savings_pct = compute_rhombic_dodecahedron_volume(box_dist, 5.0)
        st.success(
            f"💎 A rhombic dodecahedral unit cell saves ~{savings_pct:.1f}% solvent volume compared to a cubic box."
        )

    col_sim1, col_sim2, col_sim3 = st.columns(3)
    with col_sim1:
        temp_k = st.number_input("Temperature, K (Kelvin):", value=300.0, step=5.0)
    with col_sim2:
        press_bar = st.number_input("Pressure, bar (bar pressure unit):", value=1.0, step=0.1)
    with col_sim3:
        prod_ns = st.slider("Production Duration, ns (nanoseconds):", 0.01, 10.0, 0.5, 0.01)

    if prod_ns < 1.0:
        st.warning(
            "⚠️ Production duration is below 1.0 ns. Trajectories this short are "
            "generally insufficient for equilibrated stability or convergence claims; "
            "treat RMSD, RMSF, and H-bond results as exploratory only (see "
            "wiki/simulation-and-benchmark-best-practices.md)."
        )

    col_nve1, col_nve2 = st.columns([1.0, 1.0])
    with col_nve1:
        nve_enabled = st.checkbox(
            "Run NVE (Microcanonical Ensemble: Constant Particles, Volume, and Energy) drift check (LiveCoMS)",
            value=True,
            help="Short uncoupled run without thermostat or barostat to verify absence of systematic energy drift.",
        )
    with col_nve2:
        nve_ps = st.number_input(
            "NVE (Microcanonical Ensemble: Constant Particles, Volume, and Energy) Check Duration, ps (picoseconds):",
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
            "EM (Energy Minimization) structure with independent Maxwell–Boltzmann velocities and repeat "
            "NVT (Canonical Ensemble: Constant Particles, Volume, and Temperature), "
            "NPT (Isothermal–Isobaric Ensemble: Constant Particles, Pressure, and Temperature), and production."
        ),
    )
    allow_cofactor_removal = st.checkbox(
        "Allow removal of metallic cofactors (HEM/FE/MG/ZN/MN/CU/CA)",
        value=False,
        help=(
            "By default, MD (Molecular Dynamics) is refused if essential cofactors are present because "
            "removing them alters system chemistry. Enable this only after confirming the cofactor is non-essential."
        ),
    )
    retain_monatomic_ions = st.checkbox(
        "Retain monatomic metal ions (ZN/MG/CA/K)",
        value=False,
        help=(
            "When enabled, supported monatomic ions are added to complex.gro / topol.top instead of being stripped. "
            "Polyatomic cofactors such as HEM remain rejected."
        ),
    )
    retain_crystal_waters = st.checkbox(
        "Retain crystallographic water molecules (HOH/WAT)",
        value=False,
        help=(
            "Crystallographic waters are appended as SOL residues to complex.gro / topol.top. "
            "Missing hydrogens are placed geometrically using TIP3P (Transferable Intermolecular Potential with 3 Points Water Model) geometry with approximate orientation."
        ),
    )

    st.markdown(
        "**Thermostat:** `v-rescale` (Bussi Stochastic Velocity-Rescaling Thermostat) | "
        "**Coupling Groups (`tc-grps`):** `Protein_Ligand Water_and_ions`"
    )
    st.markdown(
        "**NPT (Isothermal–Isobaric Ensemble: Constant Particles, Pressure, and Temperature) Equilibration:** "
        "`berendsen` — 100 ps (picoseconds) | "
        "**Production Barostat:** `Parrinello-Rahman` | "
        "**Electrostatics:** PME (Particle Mesh Ewald Electrostatics) | "
        "**Bond Constraints:** LINCS (Linear Constraint Solver for Molecular Bonds)"
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
        "🚀 Run MD (Molecular Dynamics) Pipeline",
        type="primary",
        width="stretch",
    )
    export_clicked = col_btn2.button(
        "📥 Export Complete GROMACS (GROningen MAchine for Chemical Simulations) Protocol Bundle",
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
                with st.spinner("Building real topology via pdb2gmx and ACPYPE (AnteChamber PYthon Parser interfacE)..."):
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
                    label="⬇️ Download GROMACS (GROningen MAchine for Chemical Simulations) Bundle (.ZIP)",
                    data=zip_bytes,
                    file_name=f"md_bundle_{pose_meta.get('target_pdb', 'target')}_{ligand_resname}.zip",
                    mime="application/zip",
                    width="stretch",
                )
                st.success(
                    "Real protocol bundle ready: `topol.top` — TOP (Master System Topology File), "
                    "`complex.gro` — GRO (GROMACS Coordinate File), ligand ITP (Include Topology File) / posre, "
                    "MDP (Molecular Dynamics Parameters File) configs, `make_index.py` for NDX (Atom Index File), "
                    "and `run_md.sh`/`run_md.bat`."
                )
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
            with st.spinner("Building real topology (pdb2gmx + ACPYPE — AnteChamber PYthon Parser interfacE)..."):
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
            st.success(f"✅ Real GROMACS (GROningen MAchine for Chemical Simulations) simulation completed successfully! ({run_id})")
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
            "Mean Ligand RMSD (Root-Mean-Square Deviation)",
            f"{mean_rmsd:.2f} Å" if mean_rmsd is not None else "Not available",
            ("Stable (≤ 2.0 Å)" if mean_rmsd <= 2.0 else "Notable Fluctuation")
            if mean_rmsd is not None else None,
        )
    with col_kpi2:
        st.metric(
            "Mean H-bond (Hydrogen Bond) Count",
            f"{mean_hbonds:.1f}" if mean_hbonds is not None else "Not available",
        )
    with col_kpi3:
        production_ns = manifest.get("parameters", {}).get("production_duration_ns")
        st.metric(
            "Production Duration — ns (nanoseconds)",
            f"{production_ns:.2f} ns" if production_ns is not None else "Not recorded",
        )
    with col_kpi4:
        st.metric(
            "Engine — GROMACS (GROningen MAchine for Chemical Simulations)",
            f"GROMACS {manifest.get('engine', {}).get('version', 'unknown')}",
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
            "Within-trajectory uncertainty (RMSD — Root-Mean-Square Deviation in Å [Angstroms]; "
            "H-bond — Hydrogen Bond): " + "; ".join(uncertainty) + ". "
            f"The first {burn_in_fraction * 100:.0f}% was discarded as burn-in. "
            "These are single-trajectory error bars, not inter-replicate variance."
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
            st.caption(
                "**Across Independent Replicates:** " + " | ".join(comparisons) + ". "
                "Replicates are independent across NVT (Canonical Ensemble: Constant Particles, Volume, and Temperature), "
                "NPT (Isothermal–Isobaric Ensemble: Constant Particles, Pressure, and Temperature), and production, "
                "sharing only the EM (Energy Minimization) structure."
            )

    tab_rmsd, tab_rmsf, tab_hbond, tab_mmpbsa, tab_traj3d = st.tabs([
        "📈 Pose Stability — RMSD (Root-Mean-Square Deviation)",
        "🧬 Protein Flexibility — RMSF (Root-Mean-Square Fluctuation)",
        "🔗 Hydrogen Bonds — H-bond (Hydrogen Bond)",
        "⚡ Binding Energy — MM-PBSA (Molecular Mechanics / Poisson–Boltzmann Surface Area) / MM-GBSA (Molecular Mechanics / Generalized Born Surface Area)",
        "🎥 Trajectory Viewer (3D)",
    ])

    with tab_rmsd:
        st.plotly_chart(create_rmsd_plotly(rmsd_records, lig_name), width="stretch")
        st.markdown(r"""
        **Interpretation — RMSD (Root-Mean-Square Deviation):** For the reference Lemkul/JZ4 system the ligand heavy-atom RMSD plateaus at $\sim 1.0\text{ \AA}$ ($\text{Angstroms}$), which is the documented stability reference (see `wiki/gromacs-protein-ligand-protocol.md`). A generic upper bound of $2.0\text{ \AA}$ can indicate a stable pose, but sustained drift above $3.5\text{ \AA}$ indicates ligand dissociation or a conformational shift.
        """)

    with tab_rmsf:
        st.plotly_chart(create_rmsf_plotly(rmsf_records, md_res['pose_meta'].get('target_pdb')), width="stretch")
        st.markdown(r"""
        **Interpretation — RMSF (Root-Mean-Square Fluctuation):** High peaks correspond to flexible surface loops and termini,
        whereas low RMSF values in the active site indicate a stable interaction network with the ligand.
        """)

    with tab_hbond:
        st.plotly_chart(create_hbond_count_plotly(hbond_points, lig_name), width="stretch")
        st.caption(
            "Source: `gmx hbond -num` on the production trajectory after PBC (Periodic Boundary Conditions) re-wrapping. "
            "Criteria: 0.35 nm (nanometers) donor–acceptor distance and 30° angle."
        )

    with tab_mmpbsa:
        if mmpbsa_res is None:
            st.warning(MMPBSA_NOT_AVAILABLE_MESSAGE)
            st.code(
                "gmx_MMPBSA -O -cs prod_md.tpr -ci index.ndx -cg Protein Ligand "
                "-ct prod_md_center.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat",
                language="bash",
            )
            st.caption(
                "This command matches the official `gmx_MMPBSA` interface: `-cs` complex structure — "
                "TPR (Portable Binary Run Input File), `-ci` NDX (Atom Index File), `-cg` receptor/ligand groups, "
                "`-ct` XTC (Compressed Trajectory File), and `-cp` TOP (Master System Topology File). "
                "Run manually after verifying the trajectory is PBC (Periodic Boundary Conditions) re-wrapped."
            )
        else:
            st.plotly_chart(create_mmpbsa_breakdown_plotly(mmpbsa_res), width="stretch")
            st.markdown(
                f"**Solvation Model:** {mmpbsa_res.solvation_model} | "
                f"**ΔG_bind:** {mmpbsa_res.delta_g_bind:.2f} ± {mmpbsa_res.delta_g_bind_std:.2f} kcal/mol (kilocalories per mole)"
            )
            eps_text = (
                f"{mmpbsa_res.internal_dielectric:g}"
                if mmpbsa_res.internal_dielectric_specified
                else "Not specified"
            )
            st.caption(
                f"Frames evaluated: {mmpbsa_res.frames_evaluated or 'Not specified'} | "
                f"eps_int (internal dielectric): {eps_text} | "
                f"GB (Generalized Born) model: {mmpbsa_res.igb if mmpbsa_res.igb is not None else 'Not specified'} | "
                f"-TΔS (entropy term, kcal/mol): {mmpbsa_res.entropy_term_kcal_mol if mmpbsa_res.entropy_term_kcal_mol is not None else 'Not specified'}"
            )

    with tab_traj3d:
        st.markdown("**Final Production Frame — PBC (Periodic Boundary Conditions)-wrapped & fitted:**")
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
            st.info("Trajectory PDB (Protein Data Bank) file not found.")

    with st.expander("📁 Run Artifacts & Warnings"):
        st.json(manifest.get("artifacts", {}))
        for warning in manifest.get("warnings", []):
            st.warning(warning)
        for limitation in manifest.get("limitations", []):
            st.caption(f"⚠️ {limitation}")

    st.caption(f"🛡️ **Biophysical Disclaimer**: {DISCLAIMER_MD}")
