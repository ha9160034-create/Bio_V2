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
    "تنبيه علمي بيوفيزيائي: نتائج محاكاة الديناميكا الجزيئية تعتمد على حقول القوى الكلاسيكية "
    "وتقريب الأيونات والمذيب، بينما قيم MM-PB(GB)SA لا تُعرض إلا إذا شُغّل gmx_MMPBSA فعلياً. "
    "تُستخدم النتائج لتقييم الاستقرار وترتيب المركبات، ولا تُعد بديلاً عن القياسات المعملية (ITC / SPR)."
)


def render_mol_3d(pdb_text: str, height: int = 400, ligand_resname: str = "LIG"):
    """Render a 3D interactive molecular visualization using py3Dmol."""
    viewer = py3Dmol.view(width="100%", height=height)
    viewer.addModel(pdb_text, "pdb")
    viewer.setStyle({"hetflag": False}, {"cartoon": {"color": "spectrum", "opacity": 0.85}})
    viewer.setStyle(
        {"resn": [ligand_resname, "JZ4", "UNL"]},
        {"stick": {"colorscheme": "greenCarbon", "radius": 0.25}},
    )
    viewer.zoomTo()
    st.components.v1.html(viewer._make_html(), height=height + 20)


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
    st.header("🔬 محاكاة الديناميكا الجزيئية (Molecular Dynamics - GROMACS Pipeline)")
    st.caption("بروتوكول المحاكاة الذرية الكاملة للمركب البروتيني-الدوائي وفق معايير Lemkul و LiveCoMS الموثقة في wiki/")

    with st.expander("📖 المرجع العلمي والبروتوكول المعتمد (Wiki & Biophysical Guardrails)", expanded=False):
        st.markdown("""
        - **`wiki/gromacs-protein-ligand-protocol.md`**: بروتوكول Lemkul لإعداد وتوازن ومحاكاة المركب البروتيني-الدوائي.
        - **`wiki/open-force-field-and-ligand-parameterization.md`**: توصيف الدواء الصغير عبر ACPYPE/GAFF2 وAM1-BCC.
        - **`wiki/simulation-and-benchmark-best-practices.md`**: حظر Berendsen في الإنتاج ودمج مجموعات الترموستات.
        - **`wiki/mm-pbsa-and-mm-gbsa.md`**: MM-PB(GB)SA تُحسب فقط عبر gmx_MMPBSA على مسار الإنتاج الفعلي.
        """)

    gmx_info = get_gmx_version()
    backend = detect_gromacs()
    col_stat1, col_stat2, col_stat3 = st.columns(3)
    with col_stat1:
        if gmx_info["available"]:
            where = "WSL" if gmx_info.get("backend") == "wsl" else "Native"
            st.success(f"🟢 GROMACS متاح ({gmx_info['version']}) — {where}")
        else:
            st.error("🔴 GROMACS غير متاح. لن يتم تشغيل أي محاكاة.")
    with col_stat2:
        st.metric("تسريع العتاد (GPU)", "مُفعّل" if gmx_info.get("gpu_acceleration") else "CPU / SIMD")
    with col_stat3:
        st.metric("هندسة الصندوق الافتراضية", "اثنا عشري معيني (-23% ماء)")

    if backend is None:
        st.info(
            "💡 ثبّت GROMACS على النظام أو داخل WSL Ubuntu ثم تحقق من الأمر `gmx --version`. "
            "يمكن أيضاً ضبط المتغير `GMX_BIN`."
        )
        return

    tools = check_backend_tools(backend)
    if not tools.get("acpype"):
        st.error(
            "🔴 توصيف الدواء الصغير غير متاح: ACPYPE غير مثبت في البيئة الخلفية. "
            "ثبّته داخل WSL Ubuntu (مثال: `pip install --user acpype`)؛ ACPYPE يحمل معه "
            "ثنائيات AmberTools الخاصة به."
        )
        return

    with st.expander("🧰 فحص أدوات المحاكاة (Backend Preflight)", expanded=False):
        st.markdown(
            f"- **GROMACS:** `{gmx_info['path']}` ({gmx_info.get('backend')})\n"
            f"- **ACPYPE:** {'✅' if tools.get('acpype') else '❌'} | "
            f"**Antechamber:** {'✅ (PATH أو ACPYPE المدمج)' if tools.get('antechamber') else '❌'} | "
            f"**OpenBabel:** {'✅' if tools.get('obabel') else '❌'}\n"
            f"- **حقول القوى المتاحة:** {', '.join(tools.get('force_fields', []))}"
        )

    st.markdown("---")
    st.subheader("1️⃣ اختيار جولة الإرساء واستخلاص وضعية الارتباط المثلى")

    current_docking = st.session_state.get("current_docking_result")
    saved_runs = scan_available_docking_runs()

    run_options: List[str] = []
    run_dir_map: Dict[str, Path] = {}
    if current_docking and "run_dir" in current_docking:
        cid = current_docking.get("run_id", "active_session")
        label = (
            f"الجلسة الحالية ({cid}) — {current_docking.get('ligand_name', 'Ligand')} "
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
        st.warning("⚠️ لا توجد جولات إرساء مكتملة في `docking_runs/` أو في الجلسة الحالية.")
        st.info("💡 شغّل الإرساء أولاً من تبويب '💊 محاكاة الارتباط'، ثم عد إلى هنا. لن تُنشأ بيانات تجريبية وهمية.")
        return

    selected_run_label = st.selectbox("اختر جولة الإرساء المصدرية:", run_options)
    selected_run_dir = run_dir_map[selected_run_label]

    try:
        rec_path, best_pose, pose_meta = extract_best_pose_from_docking_run(selected_run_dir)
    except Exception as ex:  # noqa: BLE001
        st.error(f"خطأ أثناء استخراج بيانات الوضعية من الجولة: {ex}")
        return

    col_meta1, col_meta2, col_meta3, col_meta4 = st.columns(4)
    with col_meta1:
        st.metric("الهدف (Target PDB)", pose_meta.get("target_pdb", "Unknown"))
    with col_meta2:
        st.metric("اسم الدواء (Ligand)", pose_meta.get("ligand_name", "Ligand"))
    with col_meta3:
        st.metric("أفضل طاقة ارتباط (Vina Score)", f"{pose_meta['best_vina_score']:.2f} kcal/mol")
    with col_meta4:
        st.metric("عدد الذرات الثقيلة", f"{pose_meta['heavy_atoms']} (إجمالي {pose_meta['total_atoms']})")

    if not pose_meta.get("ligand_sdf_path"):
        st.warning(
            "⚠️ لم يتم العثور على ملف SDF للربيطة داخل جولة الإرساء. "
            "سيحاول النظام استنتاج رتب الروابط من الوضعية، لكن إعادة الإرساء تُنتج توصيفاً أدق."
        )

    with st.expander("👁️ معاينة وضعية الارتباط المستخلصة (3D Complex)", expanded=False):
        rec_text = Path(rec_path).read_text(encoding="utf-8", errors="replace")
        combined_preview_pdb = rec_text.replace("END\n", "") + "\n" + best_pose.pdb_block
        render_mol_3d(combined_preview_pdb, height=360, ligand_resname="LIG")

    st.markdown("---")
    st.subheader("2️⃣ توصيف المركب وحقل القوى (Topology & Force Field)")

    amber_fields = [ff for ff in tools.get("force_fields", []) if ff.startswith("amber")]
    if not amber_fields:
        st.error("لا يوجد حقل قوى Amber مثبت لهذا الـ backend. ثبّت `amber99sb-ildn.ff` أو استخدم GROMACS كامل.")
        return
    default_ff = "amber99sb-ildn" if "amber99sb-ildn" in amber_fields else amber_fields[0]

    col_top1, col_top2 = st.columns(2)
    with col_top1:
        ff_choice = st.selectbox(
            "حقل القوى (Force Field):",
            amber_fields,
            index=amber_fields.index(default_ff),
            help="المسار المدعوم حالياً: Amber99SB-ILDN للبروتين + GAFF2/AM1-BCC للربيطة.",
        )
    with col_top2:
        water_label = st.selectbox(
            "نموذج الماء الصريح (Water Model):",
            ["TIP3P", "SPC/E", "TIP4P-Ew"],
        )
    water_model = {"TIP3P": "tip3p", "SPC/E": "spc", "TIP4P-Ew": "tip4pew"}[water_label]

    st.info(
        "🔒 **المسار الفعلي الوحيد المتاح:** `gmx pdb2gmx` للبروتين + ACPYPE (GAFF2 + AM1-BCC) للربيطة. "
        "خيارات OpenFF Sage وCHARMM36/CGenFF معطّلة لأن أدواتها/ملفاتها غير مثبتة، ولن تُولَّد معاملات وهمية."
    )
    st.info(
        "🔒 **قيود المواقع الثقيلة:** تُطبَّق قيود توافقية بقوة 1000 kJ/(mol·nm²) على الذرات الثقيلة للدواء "
        "أثناء NVT وNPT لمنع انزلاقه من الجيب."
    )

    st.markdown("---")
    st.subheader("3️⃣ إعداد صندوق المحاكاة ومعاملات المحاكاة")

    col_box1, col_box2 = st.columns(2)
    with col_box1:
        box_label = st.selectbox(
            "هندسة الصندوق الحجمي (Unit Cell):",
            ["dodecahedron (اثنا عشري معيني - موصى به)", "cubic (مكعب)"],
        )
    with col_box2:
        box_dist = st.slider("المسافة العازلة (nm):", 0.8, 1.5, 1.0, 0.1)

    if "dodecahedron" in box_label:
        _, _, savings_pct = compute_rhombic_dodecahedron_volume(box_dist, 5.0)
        st.success(
            f"💎 الصندوق الاثنا عشري المعيني يوفر ~{savings_pct:.1f}% من حجم الماء مقارنة بالمكعب."
        )

    col_sim1, col_sim2, col_sim3 = st.columns(3)
    with col_sim1:
        temp_k = st.number_input("درجة الحرارة (K):", value=300.0, step=5.0)
    with col_sim2:
        press_bar = st.number_input("الضغط (bar):", value=1.0, step=0.1)
    with col_sim3:
        prod_ns = st.slider("مدة محاكاة الإنتاج (ns):", 0.01, 10.0, 0.5, 0.01)

    col_nve1, col_nve2 = st.columns([1.0, 1.0])
    with col_nve1:
        nve_enabled = st.checkbox(
            "تشغيل فحص NVE drift (LiveCoMS)",
            value=True,
            help="تشغيل قصير بدون ترموستات/باروستات للتحقق من عدم وجود انحراف منهجي في الطاقة.",
        )
    with col_nve2:
        nve_ps = st.number_input(
            "مدة فحص NVE (ps):",
            min_value=1.0,
            max_value=100.0,
            value=10.0,
            step=1.0,
            disabled=not nve_enabled,
        )
    n_replicates = st.slider(
        "عدد مسارات الإنتاج المستقلة (replicates):",
        min_value=1,
        max_value=5,
        value=3,
        help=(
            "المسار الأول هو الأساسي. المسارات الإضافية تبدأ من نفس بنية EM "
            "بسرعات Maxwell-Boltzmann مستقلة وتُعيد NVT/NPT والإنتاج."
        ),
    )
    allow_cofactor_removal = st.checkbox(
        "السماح بإزالة المرافقات المعدنية (HEM/FE/MG/ZN/MN/CU/CA)",
        value=False,
        help=(
            "افتراضياً يُرفض MD إذا وُجدت مرافقات حيوية لأن حذفها يغيّر كيمياء النظام. "
            "فعّل هذا فقط بعد التأكد من أن المرافق غير مطلوب."
        ),
    )
    retain_monatomic_ions = st.checkbox(
        "الاحتفاظ بالأيونات المعدنية الأحادية (ZN/MG/CA/K)",
        value=False,
        help=(
            "عند التفعيل تُضاف الأيونات الأحادية المدعومة إلى complex.gro/topol.top "
            "بدلاً من حذفها. المرافقات متعددة الذرات مثل HEM ما زالت مرفوضة."
        ),
    )
    retain_crystal_waters = st.checkbox(
        "الاحتفاظ بالماء البلوري (HOH/WAT)",
        value=False,
        help=(
            "يُضاف الماء البلوري كجزيئات SOL إلى complex.gro/topol.top. إذا غابت "
            "الهيدروجينات تُولّد هندسياً (TIP3P) ويكون اتجاهها تقريبياً."
        ),
    )

    st.markdown("**Thermostat:** `v-rescale` | **tc-grps:** `Protein_Ligand Water_and_ions`")
    st.markdown("**NPT equilibration:** `berendsen` (100 ps) | **Production:** `Parrinello-Rahman`")

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
        "🚀 تشغيل محاكاة الديناميكا الجزيئية (Run MD Pipeline)",
        type="primary",
        width="stretch",
    )
    export_clicked = col_btn2.button(
        "📥 تصدير حزمة GROMACS الكاملة (Export Protocol Bundle)",
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
                with st.spinner("جاري بناء طوبولوجيا حقيقية عبر pdb2gmx وACPYPE..."):
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
                    label="⬇️ تحميل حزمة GROMACS (.ZIP)",
                    data=zip_bytes,
                    file_name=f"md_bundle_{pose_meta.get('target_pdb', 'target')}_{ligand_resname}.zip",
                    mime="application/zip",
                    width="stretch",
                )
                st.success(
                    "تم تجهيز حزمة حقيقية: `topol.top`, `complex.gro`, ITP/posre للربيطة، ملفات MDP، "
                    "`make_index.py`، وسكربتات `run_md.sh`/`run_md.bat`."
                )
                if topology.warnings:
                    for w in topology.warnings:
                        st.warning(w)
            except (TopologyBuildError, Exception) as exc:  # noqa: BLE001
                st.error(f"فشل بناء حزمة GROMACS: {exc}")

    if run_clicked:
        run_id, md_run_dir = create_md_run_directory(prefix="md")
        progress_bar = st.progress(0.0)
        status_text = st.empty()

        def on_progress(stage: str, frac: float):
            progress_bar.progress(min(max(frac, 0.0), 1.0))
            status_text.markdown(f"**الحالة الآن:** {stage}")

        cfg.run_id = run_id
        try:
            with st.spinner("جاري بناء الطوبولوجيا الحقيقية (pdb2gmx + ACPYPE)..."):
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
            st.success(f"✅ اكتملت المحاكاة الفعلية عبر GROMACS بنجاح! ({run_id})")
        except (TopologyBuildError, SimulationError, Exception) as exc:  # noqa: BLE001
            _write_failure_manifest(md_run_dir, str(exc))
            st.session_state.pop("current_md_result", None)
            st.error("❌ فشل تشغيل المحاكاة. لم تُولَّد أي نتائج وهمية.")
            st.exception(exc)

    # ------------------------------------------------------------------
    # Results (only real GROMACS artifacts)
    # ------------------------------------------------------------------
    md_res = st.session_state.get("current_md_result")
    if not md_res:
        return

    manifest = md_res.get("manifest", {})
    if manifest.get("status") != "completed":
        st.error("هذه الجولة لم تكتمل بنجاح؛ لا توجد نتائج للعرض.")
        return

    st.markdown("---")
    st.subheader(
        f"📊 نتائج المسار الفعلي: {md_res['pose_meta'].get('ligand_name')} ⚯ "
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
    burn_in_fraction = float(getattr(md_res.get("config"), "analysis_burn_in_fraction", 0.2))
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
            "متوسط RMSD للدواء",
            f"{mean_rmsd:.2f} Å" if mean_rmsd is not None else "غير متاح",
            ("مستقر (≤ 2.0 Å)" if mean_rmsd <= 2.0 else "تذبذب ملحوظ")
            if mean_rmsd is not None else None,
        )
    with col_kpi2:
        st.metric("متوسط الروابط الهيدروجينية", f"{mean_hbonds:.1f}" if mean_hbonds is not None else "غير متاح")
    with col_kpi3:
        st.metric("مدة الإنتاج", f"{manifest['parameters']['production_duration_ns']:.2f} ns")
    with col_kpi4:
        st.metric("المحرك", f"GROMACS {manifest['engine']['version']}")
    uncertainty = []
    for label, summary, unit in (("RMSD", rmsd_summary, "Å"), ("H-bonds", hbond_summary, "")):
        if summary["n"]:
            uncertainty.append(
                f"{label} = {summary['mean']:.2f} ± {summary['sem']:.2f} {unit} "
                f"(block SEM {summary['block_sem']:.2f}, n={int(summary['n'])})"
            )
    if uncertainty:
        st.caption(
            "عدم اليقين داخل المسار: " + "; ".join(uncertainty) + ". "
            f"تم استبعاد أول {burn_in_fraction * 100:.0f}% كـ burn-in. "
            "هذه أشرطة خطأ داخل مسار واحد، وليست تبايناً بين replicates."
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
                "**بين المسارات المستقلة:** " + " | ".join(comparisons) + ". "
                "المسارات مستقلة في NVT/NPT/الإنتاج وتشترك فقط في بنية EM."
            )

    tab_rmsd, tab_rmsf, tab_hbond, tab_mmpbsa, tab_traj3d = st.tabs([
        "📈 استقرار الوضعية (RMSD)",
        "🧬 مرونة البروتين (RMSF)",
        "🔗 الروابط الهيدروجينية (H-Bonds)",
        "⚡ طاقة الارتباط (MM-PBSA/GBSA)",
        "🎥 عرض المسار (3D)",
    ])

    with tab_rmsd:
        st.plotly_chart(create_rmsd_plotly(rmsd_records, lig_name), width="stretch")
        st.markdown(r"""
        **التفسير (RMSD):** وصول المنحنى إلى هضبة $\le 2.0\text{ \AA}$ يدل على ثبات الوضعية.
        الارتفاع المستمر فوق $3.5\text{ \AA}$ يشير إلى انفصال الدواء أو تغير هيئوي.
        """)

    with tab_rmsf:
        st.plotly_chart(create_rmsf_plotly(rmsf_records, md_res['pose_meta'].get('target_pdb')), width="stretch")
        st.markdown(r"""
        **التفسير (RMSF):** القمم العالية تمثل الحلقات السطحية والأطراف الطرفية المرنة،
        بينما تدل قيم RMSF المنخفضة في الموقع الفعّال على شبكة تجاذب مستقرة مع الدواء.
        """)

    with tab_hbond:
        st.plotly_chart(create_hbond_count_plotly(hbond_points, lig_name), width="stretch")
        st.caption("المصدر: `gmx hbond -num` على مسار الإنتاج بعد إعادة لف PBC. المعايير: 0.35 nm و30°.")

    with tab_mmpbsa:
        if mmpbsa_res is None:
            st.warning(MMPBSA_NOT_AVAILABLE_MESSAGE)
            st.code(
                "gmx_MMPBSA -O -cs prod_md.tpr -ci index.ndx -cg Protein Ligand "
                "-ct prod_md_center.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat",
                language="bash",
            )
            st.caption(
                "هذا الأمر يطابق واجهة gmx_MMPBSA الرسمية: `-cs` بنية المركب، `-ci` الفهرس، "
                "`-cg` مجموعتا المستقبل/الربيطة، `-ct` المسار، و`-cp` الطوبولوجيا. "
                "شغّله يدوياً بعد التأكد من أن المسار مُعاد لفّه PBC."
            )
        else:
            st.plotly_chart(create_mmpbsa_breakdown_plotly(mmpbsa_res), width="stretch")
            st.markdown(
                f"**النموذج:** {mmpbsa_res.solvation_model} | "
                f"**ΔG_bind:** {mmpbsa_res.delta_g_bind:.2f} ± {mmpbsa_res.delta_g_bind_std:.2f} kcal/mol"
            )
            eps_text = (
                f"{mmpbsa_res.internal_dielectric:g}"
                if mmpbsa_res.internal_dielectric_specified
                else "غير مذكور"
            )
            st.caption(
                f"الإطارات المحللة: {mmpbsa_res.frames_evaluated or 'غير مذكورة'} | "
                f"eps_int: {eps_text} | "
                f"نموذج GB: {mmpbsa_res.igb if mmpbsa_res.igb is not None else 'غير مذكور'} | "
                f"-TΔS: {mmpbsa_res.entropy_term_kcal_mol if mmpbsa_res.entropy_term_kcal_mol is not None else 'غير مذكور'}"
            )

    with tab_traj3d:
        st.markdown("**آخر إطار من مسار الإنتاج (PBC-wrapped & fitted):**")
        last_frame = load_last_frame_pdb(run_dir / "trajectory_pbc_fit.pdb")
        if last_frame:
            render_mol_3d(last_frame, height=450, ligand_resname=md_res.get("ligand_resname", "LIG"))
        else:
            st.info("ملف المسار PDB غير موجود.")

    with st.expander("📁 ملفات الجولة والتحذيرات"):
        st.json(manifest.get("artifacts", {}))
        for warning in manifest.get("warnings", []):
            st.warning(warning)
        for limitation in manifest.get("limitations", []):
            st.caption(f"⚠️ {limitation}")

    st.caption(f"🛡️ **إخلاء مسؤولية بيوفيزيائي**: {DISCLAIMER_MD}")
