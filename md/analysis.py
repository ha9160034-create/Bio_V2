"""
Post-MD Trajectory Processing & Biophysical Analysis Layer
=========================================================

Scientific Traceability & Wiki Documentation:
---------------------------------------------
- `wiki/gromacs-protein-ligand-protocol.md`: Trajectory processing (trjconv PBC re-wrapping,
  backbone fitting), ligand heavy-atom RMSD, H-bond geometry.
- `wiki/simulation-and-benchmark-best-practices.md`: LiveCoMS Foundations in Molecular Simulations
  (Sampling convergence, burn-in trimming, statistical metrics).
- `wiki/mm-pbsa-and-mm-gbsa.md`: MM-PBSA and MM-GBSA endpoint free energy rescoring,
  Single-Trajectory approach, internal dielectric selection, and the entropy omission rule.
- `wiki/water-and-solvation.md`: Solvent-accessible surface area (SASA) non-polar models.

Biophysical Guardrails (per AGENTS.md):
---------------------------------------
1. Always apply periodic boundary condition (PBC) re-wrapping (`trjconv -pbc mol -ur compact`)
   and least-squares rotational/translational fitting onto protein backbone before calculating RMSD/RMSF.
2. Calculate ligand RMSD specifically on heavy atoms relative to the energy-minimized starting structure.
3. Hydrogen-bond criteria: Donor-acceptor distance <= 0.35 nm (3.5 Å) and angle <= 30 deg.
4. MM-PB(GB)SA calculations employ the single-trajectory approach with the -T*dS omission rule
   for congeneric series ranking to minimize numerical noise.
"""

from __future__ import annotations

import os
import re
import math
import json
import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

logger = logging.getLogger("md.analysis")


@dataclass
class RMSDRecord:
    time_ns: float
    rmsd_angstrom: float


@dataclass
class RMSFRecord:
    resnum: int
    resname: str
    rmsf_angstrom: float
    is_pocket_residue: bool = False


@dataclass
class HBondPersistence:
    donor: str
    acceptor: str
    persistence_percent: float
    avg_distance_angstrom: float
    sd_distance_angstrom: float
    status: str             # "Anchoring" | "Conserved" | "Transient"


@dataclass
class MMPBSAResult:
    delta_g_bind: float             # kcal/mol
    delta_g_bind_std: float         # kcal/mol
    delta_vdw: float                # kcal/mol
    delta_elec: float               # kcal/mol
    delta_polar_solv: float         # kcal/mol (PB or GB)
    delta_nonpolar_solv: float      # kcal/mol (SASA)
    internal_dielectric: float      # parsed eps_int; 0.0 means not specified
    solvation_model: str            # parsed model, never assumed from EGB alone
    single_trajectory: bool = True
    entropy_omitted: bool = True
    frames_evaluated: int = 0
    igb: Optional[int] = None
    entropy_term_kcal_mol: Optional[float] = None
    internal_dielectric_specified: bool = False
    source_file: str = ""


def parse_xvg_file(filepath: Path) -> List[Tuple[float, float]]:
    """Parse standard GROMACS XVG output (skipping comments and headers)."""
    filepath = Path(filepath)
    if not filepath.is_file():
        return []

    data = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith(("#", "@")):
                continue
            parts = line_str.split()
            if len(parts) >= 2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    data.append((x, y))
                except ValueError:
                    continue
    return data


def kabsch_fit(p_coords: np.ndarray, q_coords: np.ndarray) -> np.ndarray:
    """
    Superimpose coordinates p onto reference coordinates q using the Kabsch algorithm.
    Both arrays must be shaped (N, 3).
    Returns optimally rotated and translated coordinates of p.
    """
    p_cent = p_coords - p_coords.mean(axis=0)
    q_cent = q_coords - q_coords.mean(axis=0)

    # Covariance matrix
    c = np.dot(p_cent.T, q_cent)
    v, s, w = np.linalg.svd(c)

    # Reflection check
    d = (np.linalg.det(v) * np.linalg.det(w)) < 0.0
    if d:
        s[-1] = -s[-1]
        v[:, -1] = -v[:, -1]

    # Optimal rotation matrix
    r = np.dot(v, w)
    fitted = np.dot(p_cent, r) + q_coords.mean(axis=0)
    return fitted


def calculate_rmsd_from_trajectory(
    traj_pdb_path: Path,
    reference_frame_idx: int = 1,
    ligand_resnames: Tuple[str, ...] = ("LIG", "JZ4", "UNL")
) -> List[RMSDRecord]:
    """
    Parse a multi-model trajectory PDB, align the protein backbone onto the reference structure,
    and calculate ligand heavy-atom RMSD (Å) across all frames.
    """
    traj_pdb_path = Path(traj_pdb_path)
    if not traj_pdb_path.is_file():
        return []

    with open(traj_pdb_path, "r", encoding="utf-8") as f:
        content = f.read()

    models = content.split("MODEL ")
    if len(models) <= 1:
        return []

    frames_data = []
    for m in models[1:]:
        lines = m.splitlines()
        time_ns = 0.0
        for l in lines:
            if "REMARK 250 TIME:" in l:
                try:
                    time_ns = float(l.split()[3])
                except (IndexError, ValueError):
                    pass

        bb_coords = []
        lig_coords = []

        for l in lines:
            if l.startswith(("ATOM  ", "HETATM")):
                name = l[12:16].strip()
                rname = l[17:20].strip()
                elem = l[76:78].strip() if len(l) >= 78 else (name[0] if name else "C")
                is_heavy = elem.upper() != "H"
                try:
                    x = float(l[30:38].strip())
                    y = float(l[38:46].strip())
                    z = float(l[46:54].strip())
                except ValueError:
                    continue

                if name in ("CA", "N", "C") and rname not in ligand_resnames:
                    bb_coords.append((x, y, z))
                elif (rname in ligand_resnames or l.startswith("HETATM")) and is_heavy:
                    lig_coords.append((x, y, z))

        frames_data.append({
            "time_ns": time_ns,
            "bb": np.array(bb_coords),
            "lig": np.array(lig_coords)
        })

    if not frames_data:
        return []

    # Reference frame (usually frame 0 or energy minimized starting structure)
    ref_idx = min(reference_frame_idx - 1, len(frames_data) - 1)
    ref_bb = frames_data[ref_idx]["bb"]
    ref_lig = frames_data[ref_idx]["lig"]

    results: List[RMSDRecord] = []
    for idx, f in enumerate(frames_data, start=1):
        t = f["time_ns"] if f["time_ns"] > 0 else idx * 0.02
        if len(f["bb"]) == len(ref_bb) and len(f["lig"]) == len(ref_lig) and len(f["lig"]) > 0:
            # Superimpose backbone
            p_cent = f["bb"] - f["bb"].mean(axis=0)
            q_cent = ref_bb - ref_bb.mean(axis=0)
            c = np.dot(p_cent.T, q_cent)
            v, s, w = np.linalg.svd(c)
            if (np.linalg.det(v) * np.linalg.det(w)) < 0.0:
                v[:, -1] = -v[:, -1]
            r = np.dot(v, w)

            # Apply same transformation to ligand
            lig_aligned = np.dot(f["lig"] - f["bb"].mean(axis=0), r) + ref_bb.mean(axis=0)
            diff = lig_aligned - ref_lig
            rmsd_val = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
        else:
            raise ValueError(
                f"Trajectory frame {idx} has incompatible atom counts "
                f"(backbone {len(f['bb'])} vs {len(ref_bb)}, ligand {len(f['lig'])} vs {len(ref_lig)})."
            )

        results.append(RMSDRecord(time_ns=t, rmsd_angstrom=rmsd_val))

    return results


def calculate_rmsf_from_trajectory(
    traj_pdb_path: Path,
    ligand_resnames: Tuple[str, ...] = ("LIG", "JZ4", "UNL")
) -> List[RMSFRecord]:
    """
    Calculate per-residue C-alpha Root Mean Square Fluctuation (RMSF in Å)
    across all frames of a PBC-corrected trajectory.
    """
    traj_pdb_path = Path(traj_pdb_path)
    if not traj_pdb_path.is_file():
        return []

    with open(traj_pdb_path, "r", encoding="utf-8") as f:
        content = f.read()

    models = content.split("MODEL ")
    if len(models) <= 1:
        return []

    # Map residue number -> list of (x, y, z) coords across frames
    res_coords: Dict[int, List[Tuple[float, float, float]]] = {}
    res_names: Dict[int, str] = {}

    for m in models[1:]:
        for l in m.splitlines():
            if l.startswith("ATOM  "):
                name = l[12:16].strip()
                rname = l[17:20].strip()
                if name == "CA" and rname not in ligand_resnames:
                    try:
                        rnum = int(l[22:26].strip())
                        x = float(l[30:38].strip())
                        y = float(l[38:46].strip())
                        z = float(l[46:54].strip())
                        if rnum not in res_coords:
                            res_coords[rnum] = []
                            res_names[rnum] = rname
                        res_coords[rnum].append((x, y, z))
                    except ValueError:
                        continue

    records: List[RMSFRecord] = []
    for rnum in sorted(res_coords.keys()):
        arr = np.array(res_coords[rnum])
        if len(arr) > 0:
            mean_coord = arr.mean(axis=0)
            fluct = arr - mean_coord
            rmsf_val = float(np.sqrt(np.mean(np.sum(fluct ** 2, axis=1))))
        else:
            rmsf_val = 0.0

        records.append(RMSFRecord(
            resnum=rnum,
            resname=res_names.get(rnum, "ALA"),
            rmsf_angstrom=rmsf_val,
            is_pocket_residue=False
        ))

    return records


MMPBSA_NOT_AVAILABLE_MESSAGE = (
    "لم يتم العثور على نتائج MM-PB(GB)SA حقيقية. هذه القيمة لا تُحسب تقديرياً من درجة Vina؛ "
    "شغّل gmx_MMPBSA على مسار الإنتاج ثم ضع FINAL_RESULTS_MMPBSA.dat في مجلد الجولة."
)


def _first_float(text: str) -> Optional[float]:
    match = re.search(r"[-+]?\d+(?:\.\d+)?", text)
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def load_mmpbsa_results(run_dir: Path) -> Optional[MMPBSAResult]:
    """
    Parse real gmx_MMPBSA output from ``FINAL_RESULTS_MMPBSA.dat``.

    Modern gmx_MMPBSA writes Greek-delta keys (``ΔVDWAALS``, ``ΔEEL``,
    ``ΔEGB``/``ΔEPB``, ``ΔESURF``, ``ΔTOTAL``) and a ``ΔG binding`` summary.
    Earlier versions and AMBER's MMPBSA.py write the same keys without the
    prefix.  Both forms are accepted.  Only the ``Delta`` section is used for
    the binding components; otherwise the complex/receptor/ligand tables would
    be mixed together.

    ``internal_dielectric`` and the GB model are only reported when they are
    explicitly present in the output; they are never guessed from the mere
    presence of ``EGB``.

    No value is ever estimated from docking scores.
    """
    run_dir = Path(run_dir)
    candidates = [
        run_dir / "FINAL_RESULTS_MMPBSA.dat",
        run_dir / "mmpbsa" / "FINAL_RESULTS_MMPBSA.dat",
        run_dir / "FINAL_RESULTS_MMPBSA.csv",
    ]
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        return None

    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    def _normalise_key(raw: str) -> str:
        key = raw.strip().upper().replace("Δ", "").replace("δ", "")
        key = key.replace("Δ", "").strip(":")
        if key.startswith("DELTA"):
            key = key[5:].strip()
        return key.replace(" ", "")

    def _numeric_after_key(parts: list[str]) -> Optional[float]:
        for token in parts[1:]:
            candidate = token.replace(",", "").strip()
            try:
                return float(candidate)
            except ValueError:
                continue
        return None

    # Prefer the explicit delta table; fall back to the last occurrence of each
    # key for legacy files that do not print a Delta header.
    in_delta_section = False
    saw_delta_section = False
    delta_components: Dict[str, float] = {}
    legacy_components: Dict[str, float] = {}

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        lowered = line.lower()
        if "delta (complex" in lowered or lowered.startswith("deltas:") or lowered.startswith("delta:"):
            in_delta_section = True
            saw_delta_section = True
            continue

        parts = line.replace(",", " ").split()
        if not parts:
            continue
        key = _normalise_key(parts[0])
        if key in {"VDWAALS", "EEL", "EGB", "EPB", "ESURF", "TOTAL", "GGAS", "GSOLV"}:
            value = _numeric_after_key(parts)
            if value is None:
                continue
            if in_delta_section:
                delta_components[key] = value
            else:
                legacy_components[key] = value
            # A blank/separator line ends the delta energy table.
            if in_delta_section and key == "TOTAL":
                in_delta_section = False

    components = delta_components if saw_delta_section and delta_components else legacy_components
    vdw = components.get("VDWAALS")
    eel = components.get("EEL")
    egb = components.get("EGB")
    epb = components.get("EPB")
    esurf = components.get("ESURF")
    total = components.get("TOTAL")
    polar = epb if epb is not None else egb

    delta_g = None
    delta_g_std = None
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if "g binding" not in stripped.lower() and "g_bind" not in stripped.lower():
            continue
        nums = re.findall(r"[-+]?\d+(?:\.\d+)?", stripped)
        if nums:
            delta_g = float(nums[0])
            if len(nums) > 1:
                delta_g_std = float(nums[1])
            break

    # Older files may only expose ΔTOTAL; do not invent a GB decomposition.
    if delta_g is None and total is not None:
        delta_g = total

    if delta_g is None:
        total_terms = [v for v in (vdw, eel, polar, esurf) if v is not None]
        if not total_terms:
            return None
        delta_g = sum(total_terms)

    if delta_g_std is None:
        delta_g_std = 0.0

    igb_match = re.search(r"\bigb\s*=\s*(\d+)", text, flags=re.IGNORECASE)
    igb = int(igb_match.group(1)) if igb_match else None
    dielectric_match = re.search(
        r"\b(?:intdiel|eps_int|internal[_ ]dielectric)\s*[=:]\s*([-+]?\d+(?:\.\d+)?)",
        text,
        flags=re.IGNORECASE,
    )
    internal_dielectric = float(dielectric_match.group(1)) if dielectric_match else 0.0
    internal_dielectric_specified = dielectric_match is not None

    if epb is not None:
        solvation_model = "MM-PBSA (Poisson-Boltzmann)"
    elif igb is not None:
        igb_labels = {
            1: "HCT",
            2: "OBC1",
            5: "OBC2",
            7: "GBn",
            8: "GBn2",
        }
        solvation_model = f"MM-GBSA ({igb_labels.get(igb, f'igb={igb}')})"
    else:
        solvation_model = "MM-GBSA (GB model not specified in output)"

    frames_match = re.search(
        r"Calculations performed using\s+(\d+)\s+complex frames",
        text,
        flags=re.IGNORECASE,
    )
    frames_evaluated = int(frames_match.group(1)) if frames_match else 0

    entropy_term = None
    if total is not None:
        entropy_term = round(float(delta_g) - float(total), 2)
    entropy_omitted = entropy_term is None or abs(entropy_term) < 1e-6
    single_trajectory = "multiple trajectory" not in text.lower()

    return MMPBSAResult(
        delta_g_bind=round(float(delta_g), 2),
        delta_g_bind_std=round(float(delta_g_std), 2),
        delta_vdw=round(float(vdw), 2) if vdw is not None else 0.0,
        delta_elec=round(float(eel), 2) if eel is not None else 0.0,
        delta_polar_solv=round(float(polar), 2) if polar is not None else 0.0,
        delta_nonpolar_solv=round(float(esurf), 2) if esurf is not None else 0.0,
        internal_dielectric=internal_dielectric,
        solvation_model=solvation_model,
        single_trajectory=single_trajectory,
        entropy_omitted=entropy_omitted,
        frames_evaluated=frames_evaluated,
        igb=igb,
        entropy_term_kcal_mol=entropy_term,
        internal_dielectric_specified=internal_dielectric_specified,
        source_file=path.name,
    )


def summarize_timeseries(
    points: List[Tuple[float, float]],
    burn_in_fraction: float = 0.2,
    n_blocks: int = 5,
) -> Dict[str, float]:
    """Burn-in-aware mean and uncertainty for a time series.

    Returns the arithmetic mean, sample standard deviation, standard error of
    the mean, and a block-averaged SEM (preferred when frames are correlated).
    A single trajectory still cannot provide replicate-to-replicate variance;
    this only quantifies within-trajectory uncertainty.
    """
    if not points:
        return {"n": 0.0, "mean": 0.0, "sd": 0.0, "sem": 0.0, "block_sem": 0.0}
    ordered = sorted(points, key=lambda point: point[0])
    latest = ordered[-1][0]
    fraction = min(max(float(burn_in_fraction), 0.0), 0.95)
    if latest > 0 and fraction > 0:
        ordered = [point for point in ordered if point[0] >= fraction * latest]
    if not ordered:
        ordered = [points[-1]]
    values = [float(point[1]) for point in ordered]
    n = len(values)
    mean = sum(values) / n
    sd = math.sqrt(sum((value - mean) ** 2 for value in values) / (n - 1)) if n > 1 else 0.0
    sem = sd / math.sqrt(n) if n > 1 else 0.0
    block_sem = 0.0
    if n >= n_blocks >= 2:
        block_size = n // n_blocks
        block_means = [
            sum(values[i * block_size:(i + 1) * block_size]) / block_size
            for i in range(n_blocks)
        ]
        block_mean = sum(block_means) / n_blocks
        block_sd = math.sqrt(
            sum((value - block_mean) ** 2 for value in block_means) / (n_blocks - 1)
        )
        block_sem = block_sd / math.sqrt(n_blocks)
    return {
        "n": float(n),
        "mean": mean,
        "sd": sd,
        "sem": sem,
        "block_sem": block_sem,
    }


def create_rmsd_plotly(rmsd_records: List[RMSDRecord], ligand_name: str = "Ligand") -> go.Figure:
    """Generate interactive Plotly visualization for Ligand Heavy-Atom RMSD."""
    if not rmsd_records:
        fig = go.Figure()
        fig.add_annotation(text="لا توجد بيانات RMSD متاحة", showarrow=False)
        return fig

    times = [r.time_ns for r in rmsd_records]
    vals = [r.rmsd_angstrom for r in rmsd_records]
    mean_val = float(np.mean(vals))

    fig = go.Figure()
    # Main RMSD curve
    fig.add_trace(go.Scatter(
        x=times,
        y=vals,
        mode="lines",
        name=f"RMSD: {ligand_name}",
        line=dict(color="#2563EB", width=2.5),
        hovertemplate="Time: %{x:.3f} ns<br>RMSD: %{y:.2f} Å<extra></extra>"
    ))

    # Running average or mean line
    fig.add_trace(go.Scatter(
        x=[times[0], times[-1]],
        y=[mean_val, mean_val],
        mode="lines",
        name=f"المتوسط ({mean_val:.2f} Å)",
        line=dict(color="#DC2626", width=1.5, dash="dash")
    ))

    # Stability threshold zone (<= 2.0 Å is considered stable pose per Lemkul protocol)
    fig.add_hrect(
        y0=0.0, y1=2.0,
        fillcolor="rgba(34, 197, 94, 0.12)",
        line_width=0,
        annotation_text="نطاق الاستقرار الهيكلي (RMSD <= 2.0 Å)",
        annotation_position="bottom right"
    )

    fig.update_layout(
        title=f"<b>الانحراف المعياري للذرات الثقيلة (Heavy-Atom RMSD)</b> — {ligand_name}",
        xaxis_title="زمن المحاكاة (Simulation Time, ns)",
        yaxis_title="RMSD (Å)",
        template="plotly_white",
        hovermode="x unified",
        margin=dict(l=40, r=40, t=60, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    )
    return fig


def create_rmsf_plotly(rmsf_records: List[RMSFRecord], target_id: str = "Target") -> go.Figure:
    """Generate interactive Plotly visualization for per-residue protein RMSF."""
    if not rmsf_records:
        fig = go.Figure()
        fig.add_annotation(text="لا توجد بيانات RMSF متاحة", showarrow=False)
        return fig

    resnums = [r.resnum for r in rmsf_records]
    vals = [r.rmsf_angstrom for r in rmsf_records]
    names = [f"{r.resname}{r.resnum}" for r in rmsf_records]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=resnums,
        y=vals,
        name="RMSF (Cα)",
        marker=dict(
            color=vals,
            colorscale="Viridis",
            showscale=True,
            colorbar=dict(title="Å")
        ),
        text=names,
        hovertemplate="Residue: %{text}<br>RMSF: %{y:.2f} Å<extra></extra>"
    ))

    fig.update_layout(
        title=f"<b>تذبذب الثمالات البروتينية (Protein Residue RMSF)</b> — {target_id}",
        xaxis_title="رقم الثمالة (Residue Index)",
        yaxis_title="RMSF (Å)",
        template="plotly_white",
        margin=dict(l=40, r=40, t=60, b=40)
    )
    return fig


def create_hbond_persistence_plotly(hbonds: List[Dict[str, Any]]) -> go.Figure:
    """Generate interactive Plotly visualization for Hydrogen-Bond Persistence."""
    if not hbonds:
        fig = go.Figure()
        fig.add_annotation(text="لا توجد بيانات روابط هيدروجينية متاحة", showarrow=False)
        return fig

    labels = [f"{h['donor']} ➔ {h['acceptor']}" for h in hbonds]
    pers = [h["persistence_percent"] for h in hbonds]
    colors = ["#16A34A" if p >= 75.0 else ("#D97706" if p >= 40.0 else "#94A3B8") for p in pers]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=pers,
        y=labels,
        orientation="h",
        marker=dict(color=colors),
        hovertemplate="التفاعل: %{y}<br>الاستدامة: %{x:.1f}%<extra></extra>"
    ))

    fig.update_layout(
        title="<b>استدامة الروابط الهيدروجينية (Hydrogen-Bond Persistence %)</b>",
        xaxis_title="نسبة الثبات عبر مسار المحاكاة (%)",
        yaxis_title="الزوج المتفاعل (Donor ➔ Acceptor)",
        xaxis=dict(range=[0, 105]),
        template="plotly_white",
        margin=dict(l=40, r=40, t=60, b=40)
    )
    return fig


def create_hbond_count_plotly(
    points: List[Tuple[float, float]],
    ligand_name: str = "Ligand",
) -> go.Figure:
    """Plot the real protein–ligand H-bond count from ``gmx hbond -num``."""
    fig = go.Figure()
    if not points:
        fig.add_annotation(text="لا توجد بيانات روابط هيدروجينية (gmx hbond)", showarrow=False)
        fig.update_layout(template="plotly_white")
        return fig

    times = [p[0] for p in points]
    counts = [p[1] for p in points]
    mean_count = sum(counts) / len(counts)

    fig.add_trace(go.Scatter(
        x=times,
        y=counts,
        mode="lines",
        name=f"H-bonds: Protein ⚯ {ligand_name}",
        line=dict(color="#059669", width=2.0),
        hovertemplate="Time: %{x:.1f} ps<br>H-bonds: %{y:.0f}<extra></extra>",
    ))
    fig.add_hline(
        y=mean_count,
        line_dash="dash",
        line_color="#0F172A",
        annotation_text=f"المتوسط: {mean_count:.1f}",
    )
    fig.update_layout(
        title="<b>عدد الروابط الهيدروجينية بين البروتين والدواء عبر الزمن</b>",
        xaxis_title="الزمن (ps)",
        yaxis_title="عدد الروابط",
        template="plotly_white",
        margin=dict(l=40, r=40, t=60, b=40),
    )
    return fig


def load_last_frame_pdb(trajectory_pdb: Path) -> str:
    """Return the last MODEL of a multi-model PDB trajectory."""
    trajectory_pdb = Path(trajectory_pdb)
    if not trajectory_pdb.is_file():
        return ""
    text = trajectory_pdb.read_text(encoding="utf-8", errors="replace")
    if "MODEL" not in text:
        return text
    blocks = [b for b in text.split("MODEL") if "ATOM" in b or "HETATM" in b]
    if not blocks:
        return ""
    last = blocks[-1]
    lines = []
    for line in last.splitlines():
        if line.startswith(("ATOM", "HETATM", "TER", "END")):
            lines.append(line)
    return "\n".join(lines) + "\nEND\n"


def create_mmpbsa_breakdown_plotly(result: MMPBSAResult) -> go.Figure:
    """Generate interactive Plotly waterfall / bar chart for MM-PB(GB)SA energy decomposition."""
    terms = [
        "فان دير فالس (ΔE_vdw)",
        "الكولومبية (ΔE_elec)",
        "الذوبان القطبي (ΔG_polar)",
        "الذوبان غير القطبي (ΔG_nonpolar)",
        "صافي طاقة الارتباط (ΔG_bind)"
    ]
    vals = [
        result.delta_vdw,
        result.delta_elec,
        result.delta_polar_solv,
        result.delta_nonpolar_solv,
        result.delta_g_bind
    ]
    colors = ["#2563EB", "#7C3AED", "#DC2626", "#059669", "#0F172A"]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=terms,
        y=vals,
        marker=dict(color=colors),
        text=[f"{v:.1f} kcal/mol" for v in vals],
        textposition="auto",
        hovertemplate="%{x}: %{y:.2f} kcal/mol<extra></extra>"
    ))

    fig.update_layout(
        title=f"<b>تحليل مكونات طاقة الارتباط الحرة ({result.solvation_model})</b>",
        yaxis_title="الطاقة (kcal/mol)",
        template="plotly_white",
        margin=dict(l=40, r=40, t=60, b=40)
    )
    return fig

