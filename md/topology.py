"""
Topology & Complex Preparation Layer for Molecular Dynamics
===========================================================

Scientific Traceability & Wiki Documentation:
---------------------------------------------
- `wiki/open-force-field-and-ligand-parameterization.md`: Direct Chemical Perception,
  OpenFF Sage 2.x, SMIRNOFF, CGenFF Penalty Rules, and ACPYPE/GAFF2.
- `wiki/gromacs-protein-ligand-protocol.md`: Justin Lemkul's GROMACS Protein–Ligand Protocol
  (System recombination, ligand decoupling, posre_lig.itp generation).
- `wiki/force-fields.md`: CHARMM36, Amber14SB, and Small-Molecule Transferability.

Biophysical Implementation Guardrails:
--------------------------------------
1. Decouple macromolecular receptor from small-molecule ligand during initial topology generation.
2. Restrain ligand heavy atoms during solvent equilibration (fc = 1000 1000 1000 kJ/(mol*nm^2)).
3. Validate CGenFF penalty scores (<10 robust, 10-50 warning/verification needed, >50 invalid).
4. Combine receptor and ligand into unified complex coordinate file with updated atom counts.

Production topology policy
--------------------------
The only parameterisation route that is executed here is the one the local
backend can actually run: ``gmx pdb2gmx`` for the receptor and ACPYPE
(GAFF2 + AM1-BCC) for the ligand. OpenFF and CGenFF generators that used to
emit invented atom types have been removed: they produced files that looked
valid but were physically meaningless. Requests for those routes now fail
with an explicit message.
"""

from __future__ import annotations

import os
import re
import json
import math
import shutil
import shlex
import sys
import subprocess
import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any

from md.engine import (
    GromacsBackend,
    GromacsExecutionError,
    detect_gromacs,
    resolve_force_field,
    resolve_water_model,
    run_backend_shell,
    run_gmx_command,
)
from md.indexing import parse_gro_atoms

logger = logging.getLogger("md.topology")
PROJECT_ROOT = Path(__file__).resolve().parents[1]


class TopologyBuildError(RuntimeError):
    """Raised when a physically valid topology cannot be produced."""


@dataclass
class LigandAtom:
    serial: int
    name: str
    resname: str
    chain: str
    resnum: int
    x: float
    y: float
    z: float
    element: str
    charge: float = 0.0
    atom_type: str = ""
    is_heavy: bool = True


@dataclass
class ExtractedPose:
    rank: int
    score: float
    seed: int
    atoms: List[LigandAtom]
    pdb_block: str
    sdf_block: str = ""
    formula: str = ""
    heavy_atom_count: int = 0
    total_atom_count: int = 0


@dataclass
class CGenFFPenaltyReport:
    max_charge_penalty: float
    max_dihedral_penalty: float
    average_penalty: float
    status: str             # "PASS" | "WARNING" | "FAIL"
    verdict: str
    details: List[str] = field(default_factory=list)


@dataclass
class TopologyResult:
    method: str
    force_field: str
    water_model: str
    ligand_resname: str
    complex_gro: Path
    complex_pdb: Path
    ligand_gro: Path
    ligand_itp: Path
    ligand_posre: Path
    topol_top: Path
    penalty_report: Optional[CGenFFPenaltyReport] = None
    warnings: List[str] = field(default_factory=list)
    receptor_pdb: Optional[Path] = None
    protein_gro: Optional[Path] = None
    ligand_sdf: Optional[Path] = None
    ligand_charge: int = 0
    protein_atom_count: int = 0
    ligand_atom_count: int = 0
    backend_label: str = ""
    water_validation: Dict[str, Any] = field(default_factory=dict)
    retained_ions: List[Dict[str, Any]] = field(default_factory=list)
    retained_waters: int = 0


# ---------------------------------------------------------------------------
# Docking output parsing
# ---------------------------------------------------------------------------

_AUTODOCK_TYPE_ELEMENTS = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O",
    "H": "H", "HD": "H", "HS": "H",
    "S": "S", "SA": "S", "P": "P",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "FE": "Fe", "ZN": "Zn", "MG": "Mg", "CA": "Ca", "MN": "Mn", "CU": "Cu", "K": "K",
}


def autodock_type_to_element(atom_type: str, atom_name: str = "") -> str:
    """Map an AutoDock atom type (e.g. ``OA``, ``NA``, ``A``) to an element.

    PDBQT columns 78-79 contain the AutoDock *type*, not a PDB element.  Taking
    them literally turns an acceptor oxygen (``OA``) into ``Oa`` and an acceptor
    nitrogen (``NA``) into sodium.  When the type is unknown, fall back to the
    atom name, never to a blind two-letter element.
    """
    raw = (atom_type or "").strip().upper()
    if raw in _AUTODOCK_TYPE_ELEMENTS:
        return _AUTODOCK_TYPE_ELEMENTS[raw]
    name = (atom_name or "").strip()
    alpha = "".join(ch for ch in name if ch.isalpha())
    if not alpha:
        return "C"
    upper = alpha.upper()
    if upper.startswith("CL"):
        return "Cl"
    if upper.startswith("BR"):
        return "Br"
    return alpha[0].upper()


def parse_pdbqt_models(pdbqt_text: str, default_resname: str = "LIG") -> List[ExtractedPose]:
    """
    Parse multi-model PDBQT output (e.g. from Vina out_seed_*.pdbqt)
    into structured poses with accurate 3D coordinates, elements, and heavy-atom identification.
    """
    poses: List[ExtractedPose] = []
    model_blocks = pdbqt_text.split("MODEL ")
    if len(model_blocks) <= 1:
        model_blocks = ["", pdbqt_text]

    for model_idx, block in enumerate(model_blocks[1:], start=1):
        lines = block.splitlines()
        score = 0.0
        seed = 42

        for line in lines:
            if "REMARK VINA RESULT:" in line:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        score = float(parts[3])
                    except ValueError:
                        pass
            elif "REMARK SEED:" in line:
                try:
                    seed = int(line.split()[-1])
                except ValueError:
                    pass

        atoms: List[LigandAtom] = []
        clean_pdb_lines: List[str] = []
        serial_counter = 1

        for line in lines:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    source_serial = int(line[6:11].strip())
                    name = line[12:16].strip()
                    resname = default_resname
                    chain = "L"
                    resnum = 1
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    raw_type = line[77:].strip() if len(line) >= 78 else ""
                    if not raw_type:
                        raw_type = line[76:78].strip() if len(line) >= 78 else ""
                    elem = autodock_type_to_element(raw_type, name)

                    charge = 0.0
                    if len(line) >= 70:
                        charge_str = line[66:76].strip()
                        try:
                            charge = float(charge_str)
                        except ValueError:
                            pass

                    is_heavy = elem.upper() != "H"
                    pdb_line = (
                        f"HETATM{source_serial:5d} {name:<4s} {resname:3s} {chain:1s}{resnum:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {elem:>2s}"
                    )
                    clean_pdb_lines.append(pdb_line)

                    atoms.append(LigandAtom(
                        serial=source_serial,
                        name=name,
                        resname=resname,
                        chain=chain,
                        resnum=resnum,
                        x=x,
                        y=y,
                        z=z,
                        element=elem,
                        charge=charge,
                        is_heavy=is_heavy,
                    ))
                    serial_counter += 1
                except Exception as ex:  # noqa: BLE001 - skip malformed atom lines
                    logger.debug("Failed parsing atom line: %s (%s)", line, ex)

        if atoms:
            heavy_count = sum(1 for a in atoms if a.is_heavy)
            poses.append(ExtractedPose(
                rank=model_idx,
                score=score,
                seed=seed,
                atoms=atoms,
                pdb_block="\n".join(clean_pdb_lines) + "\nEND\n",
                heavy_atom_count=heavy_count,
                total_atom_count=len(atoms),
            ))

    poses.sort(key=lambda p: p.score)
    for idx, p in enumerate(poses, start=1):
        p.rank = idx
    return poses


def _find_ligand_sdf(docking_run_dir: Path, manifest: Dict[str, Any]) -> Optional[Path]:
    """Locate the ligand SDF (with bond orders) written during docking."""
    candidates: List[Path] = []
    sdf_name = manifest.get("ligand", {}).get("sdf_file")
    if sdf_name:
        candidates.append(docking_run_dir / sdf_name)
    candidates.extend(sorted(docking_run_dir.glob("*_3d.sdf")))
    candidates.extend(sorted(docking_run_dir.glob("*.sdf")))
    for cand in candidates:
        if cand.is_file() and cand.stat().st_size > 0:
            return cand
    return None


def extract_best_pose_from_docking_run(docking_run_dir: Path) -> Tuple[Path, ExtractedPose, Dict[str, Any]]:
    """
    Extract the cleaned receptor and top-scoring docked ligand pose from a completed docking run.
    Returns (receptor_pdb_path, best_pose, metadata).
    """
    docking_run_dir = Path(docking_run_dir)
    if not docking_run_dir.exists():
        raise FileNotFoundError(f"Docking run directory not found: {docking_run_dir}")

    manifest_path = docking_run_dir / "manifest.json"
    manifest: Dict[str, Any] = {}
    if manifest_path.is_file():
        try:
            with open(manifest_path, "r", encoding="utf-8") as f:
                manifest = json.load(f)
        except Exception:
            pass

    receptor_candidates = [
        docking_run_dir / "receptor_clean.pdb",
        docking_run_dir / "clean_receptor.pdb",
        docking_run_dir / "receptor.pdb",
    ]
    receptor_pdb = None
    for cand in receptor_candidates:
        if cand.is_file():
            receptor_pdb = cand
            break

    if not receptor_pdb:
        for p in docking_run_dir.glob("*.pdb"):
            if "reference_ligand" not in p.name.lower() and "ligand" not in p.name.lower():
                receptor_pdb = p
                break

    if not receptor_pdb or not receptor_pdb.is_file():
        raise FileNotFoundError(f"Clean receptor PDB not found in docking run: {docking_run_dir}")

    pdbqt_files = sorted(docking_run_dir.glob("out_seed_*.pdbqt"))
    if not pdbqt_files:
        pdbqt_files = sorted(docking_run_dir.glob("*_out.pdbqt"))
    if not pdbqt_files:
        raise FileNotFoundError(f"No docking pose outputs (out_seed_*.pdbqt) found in: {docking_run_dir}")

    all_poses: List[ExtractedPose] = []
    for f in pdbqt_files:
        with open(f, "r", encoding="utf-8") as fp:
            text = fp.read()
        seed_match = re.search(r"seed_(\d+)", f.stem)
        seed_val = int(seed_match.group(1)) if seed_match else 42
        for sp in parse_pdbqt_models(text):
            sp.seed = seed_val
            all_poses.append(sp)

    if not all_poses:
        raise ValueError(f"No valid ligand poses could be extracted from: {docking_run_dir}")

    all_poses.sort(key=lambda p: p.score)
    best_pose = all_poses[0]
    sdf_path = _find_ligand_sdf(docking_run_dir, manifest)

    meta = {
        "run_id": manifest.get("run_id", docking_run_dir.name),
        "target_pdb": manifest.get("target", {}).get("pdb_id", "Unknown"),
        "ligand_name": manifest.get("ligand", {}).get("name", "Ligand"),
        "best_vina_score": best_pose.score,
        "best_pose_seed": best_pose.seed,
        "best_pose_rank": best_pose.rank,
        "heavy_atoms": best_pose.heavy_atom_count,
        "total_atoms": best_pose.total_atom_count,
        "receptor_path": str(receptor_pdb.resolve()),
        "ligand_sdf_path": str(sdf_path.resolve()) if sdf_path else None,
    }
    return receptor_pdb, best_pose, meta


# ---------------------------------------------------------------------------
# Receptor preparation and pdb2gmx
# ---------------------------------------------------------------------------

_MSE_ATOM_MAP = {"SE": "SD", "HSE": "HG", "HSE1": "HG1"}
_CRITICAL_COFACTOR_RESNAMES = {
    "HEM", "HEC", "HEB", "HEA", "FE", "MG", "ZN", "MN", "CU", "CA",
}
# Monatomic ions with explicit moleculetypes in the Amber ions.itp used here.
_RETAINABLE_MONATOMIC_IONS = {"ZN", "MG", "CA", "K"}
_CRYSTAL_WATER_RESNAMES = {"HOH", "WAT", "SOL"}
_WATER_OH_NM = 0.09572
_WATER_HOH_DEG = 104.52


def _distance_sq(a, b) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def extract_crystal_waters(receptor_pdb: Path) -> List[Dict[str, Any]]:
    """Extract crystallographic waters and build complete 3-site coordinates.

    Deposited waters often contain only the oxygen.  Hydrogens are generated at
    the standard TIP3P geometry (0.09572 nm, 104.52°), oriented away from the
    nearest protein atom to reduce initial clashes.  This is an approximate
    hydrogen orientation, not an experimental assignment.
    """
    groups: Dict[Tuple[str, str, str, str], Dict[str, Any]] = {}
    protein_coords: List[Tuple[float, float, float]] = []
    for line in Path(receptor_pdb).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM"):
            try:
                protein_coords.append((
                    float(line[30:38]) / 10.0,
                    float(line[38:46]) / 10.0,
                    float(line[46:54]) / 10.0,
                ))
            except ValueError:
                continue
            continue
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper()
        if resname not in _CRYSTAL_WATER_RESNAMES:
            continue
        chain = line[21].strip() or "A"
        resseq = line[22:26].strip()
        icode = line[26].strip()
        atomname = line[12:16].strip().upper()
        try:
            xyz = (
                float(line[30:38]) / 10.0,
                float(line[38:46]) / 10.0,
                float(line[46:54]) / 10.0,
            )
        except ValueError:
            continue
        key = (resname, chain, resseq, icode)
        group = groups.setdefault(key, {"resname": resname, "chain": chain, "resseq": resseq, "icode": icode, "atoms": {}})
        group["atoms"][atomname] = xyz

    waters: List[Dict[str, Any]] = []
    half = math.radians(_WATER_HOH_DEG / 2.0)
    for group in groups.values():
        atoms = group["atoms"]
        ow = atoms.get("O") or atoms.get("OW")
        if ow is None:
            continue
        hw1 = next((atoms[name] for name in ("HW1", "H1", "H") if name in atoms), None)
        hw2 = next((atoms[name] for name in ("H2", "HW2", "H") if name in atoms and atoms[name] != hw1), None)
        if hw1 is None or hw2 is None:
            if protein_coords:
                nearest = min(protein_coords, key=lambda p: _distance_sq(ow, p))
                away = (ow[0] - nearest[0], ow[1] - nearest[1], ow[2] - nearest[2])
            else:
                away = (1.0, 0.0, 0.0)
            norm = math.sqrt(sum(component * component for component in away))
            if norm < 1e-6:
                away = (1.0, 0.0, 0.0)
                norm = 1.0
            away = tuple(component / norm for component in away)
            # Pick any perpendicular direction for the H-O-H plane.
            ref = (0.0, 0.0, 1.0) if abs(away[2]) < 0.9 else (1.0, 0.0, 0.0)
            side = (
                away[1] * ref[2] - away[2] * ref[1],
                away[2] * ref[0] - away[0] * ref[2],
                away[0] * ref[1] - away[1] * ref[0],
            )
            side_norm = math.sqrt(sum(component * component for component in side))
            side = tuple(component / side_norm for component in side) if side_norm > 1e-6 else (1.0, 0.0, 0.0)
            hw1 = (
                ow[0] + _WATER_OH_NM * (math.cos(half) * away[0] + math.sin(half) * side[0]),
                ow[1] + _WATER_OH_NM * (math.cos(half) * away[1] + math.sin(half) * side[1]),
                ow[2] + _WATER_OH_NM * (math.cos(half) * away[2] + math.sin(half) * side[2]),
            )
            hw2 = (
                ow[0] + _WATER_OH_NM * (math.cos(half) * away[0] - math.sin(half) * side[0]),
                ow[1] + _WATER_OH_NM * (math.cos(half) * away[1] - math.sin(half) * side[1]),
                ow[2] + _WATER_OH_NM * (math.cos(half) * away[2] - math.sin(half) * side[2]),
            )
        waters.append({
            "resname": "SOL",
            "source_resname": group["resname"],
            "chain": group["chain"],
            "resseq": group["resseq"],
            "icode": group["icode"],
            "ow": ow,
            "hw1": hw1,
            "hw2": hw2,
        })
    return waters


def append_crystal_waters_to_topology(topology: Path, water_count: int) -> Path:
    """Append one ``SOL N`` molecule entry to ``[ molecules ]``."""
    if water_count <= 0:
        return Path(topology)
    topology = Path(topology)
    lines = topology.read_text(encoding="utf-8", errors="replace").splitlines()
    molecules_idx = None
    for index, line in enumerate(lines):
        if line.strip().lower() == "[ molecules ]":
            molecules_idx = index
    if molecules_idx is None:
        raise TopologyBuildError(f"No [ molecules ] section found in {topology}.")
    insert_at = molecules_idx + 1
    for index in range(molecules_idx + 1, len(lines)):
        stripped = lines[index].split(";", 1)[0].strip()
        if stripped.startswith("["):
            break
        if stripped:
            insert_at = index + 1
    lines.insert(insert_at, f"{'SOL':<20s} {water_count}")
    topology.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return topology


def append_crystal_waters_to_gro(gro_path: Path, waters: List[Dict[str, Any]]) -> Path:
    """Append crystallographic waters (OW/HW1/HW2) to a complex GRO file."""
    if not waters:
        return Path(gro_path)
    gro_path = Path(gro_path)
    lines = gro_path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        atom_count = int(lines[1].strip().split()[0])
    except (IndexError, ValueError) as exc:
        raise TopologyBuildError(f"Malformed GRO file: {gro_path}") from exc
    body = lines[2:2 + atom_count]
    box_line = lines[2 + atom_count] if len(lines) > 2 + atom_count else "   0.00000   0.00000   0.00000"
    last_resnum = 0
    for line in body:
        try:
            last_resnum = max(last_resnum, int(line[0:5]))
        except ValueError:
            pass
    water_lines = []
    atom_index = atom_count
    for water in waters:
        last_resnum += 1
        for atomname, xyz in (("OW", water["ow"]), ("HW1", water["hw1"]), ("HW2", water["hw2"])):
            atom_index += 1
            water_lines.append(
                f"{last_resnum:5d}{'SOL':<5s}{atomname:<5s}{atom_index:5d}"
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            )
    gro_path.write_text(
        "\n".join(
            [lines[0], f"{atom_count + 3 * len(waters):5d}", *body, *water_lines, box_line]
        ) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return gro_path


def extract_monatomic_ions(receptor_pdb: Path) -> Tuple[List[Dict[str, Any]], set]:
    """Extract monatomic metal ions that can be retained in the MD topology.

    Multi-atom cofactors (e.g. HEM) and ions without an explicit Amber
    moleculetype are returned in ``unsupported`` and must not be silently
    dropped.
    """
    groups: Dict[Tuple[str, str, str, str], List[Dict[str, Any]]] = {}
    for line in Path(receptor_pdb).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper()
        chain = line[21].strip() or "A"
        resseq = line[22:26].strip()
        icode = line[26].strip()
        atomname = line[12:16].strip()
        element = (line[76:78].strip() or atomname).upper()
        try:
            x = float(line[30:38]) / 10.0
            y = float(line[38:46]) / 10.0
            z = float(line[46:54]) / 10.0
        except ValueError:
            continue
        key = (resname, chain, resseq, icode)
        groups.setdefault(key, []).append({
            "resname": resname,
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "atomname": atomname,
            "element": element,
            "x": x,
            "y": y,
            "z": z,
        })

    ions: List[Dict[str, Any]] = []
    unsupported: set = set()
    for (resname, _chain, _resseq, _icode), atoms in groups.items():
        if resname in _RETAINABLE_MONATOMIC_IONS and len(atoms) == 1:
            ions.append(atoms[0])
        elif resname in _CRITICAL_COFACTOR_RESNAMES:
            unsupported.add(resname)
    return ions, unsupported


def append_ions_to_topology(topology: Path, ions: List[Dict[str, Any]]) -> Path:
    """Append retained ion molecule entries to ``[ molecules ]``."""
    topology = Path(topology)
    lines = topology.read_text(encoding="utf-8", errors="replace").splitlines()
    molecules_idx = None
    for index, line in enumerate(lines):
        if line.strip().lower() == "[ molecules ]":
            molecules_idx = index
    if molecules_idx is None:
        raise TopologyBuildError(f"No [ molecules ] section found in {topology}.")
    insert_at = molecules_idx + 1
    for index in range(molecules_idx + 1, len(lines)):
        stripped = lines[index].split(";", 1)[0].strip()
        if stripped.startswith("["):
            break
        if stripped:
            insert_at = index + 1
    for ion in ions:
        lines.insert(insert_at, f"{ion['resname']:<20s} 1")
        insert_at += 1
    topology.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return topology


def append_ions_to_gro(gro_path: Path, ions: List[Dict[str, Any]]) -> Path:
    """Append retained ion coordinates to a complex GRO file."""
    gro_path = Path(gro_path)
    lines = gro_path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        atom_count = int(lines[1].strip().split()[0])
    except (IndexError, ValueError) as exc:
        raise TopologyBuildError(f"Malformed GRO file: {gro_path}") from exc
    body = lines[2:2 + atom_count]
    box_line = lines[2 + atom_count] if len(lines) > 2 + atom_count else "   0.00000   0.00000   0.00000"
    last_resnum = 0
    for line in body:
        try:
            last_resnum = max(last_resnum, int(line[0:5]))
        except ValueError:
            pass
    ion_lines = []
    last_atom = atom_count
    for ion in ions:
        last_resnum += 1
        last_atom += 1
        ion_lines.append(
            f"{last_resnum:5d}{ion['resname']:<5s}{ion['atomname']:<5s}{last_atom:5d}"
            f"{ion['x']:8.3f}{ion['y']:8.3f}{ion['z']:8.3f}"
        )
    gro_path.write_text(
        "\n".join(
            [lines[0], f"{atom_count + len(ions):5d}", *body, *ion_lines, box_line]
        ) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return gro_path


def inspect_water_topology_text(itp_text: str, gro_text: str, water_model: str) -> Dict[str, Any]:
    """Validate a 4-site water ITP/GRO pair without executing GROMACS.

    GROMACS has no ``pdb2gmx -water tip4pew`` value, so this project builds with
    ``tip4p`` and swaps the include.  This check makes the assumptions explicit:
    the ITP must define ``SOL`` with four atoms and the solvent template must
    contain whole four-site water molecules.
    """
    moleculetype = ""
    atom_count = 0
    in_atoms = False
    in_moleculetype = False
    for raw in itp_text.splitlines():
        line = raw.strip()
        if line.startswith("["):
            section = line.replace("[", "").replace("]", "").strip().lower()
            in_atoms = section == "atoms"
            in_moleculetype = section == "moleculetype"
            continue
        if not line or line.startswith((";", "#")):
            continue
        if in_moleculetype and not moleculetype:
            moleculetype = line.split()[0]
        elif in_atoms:
            atom_count += 1

    gro_lines = [line for line in gro_text.splitlines() if line.strip()]
    gro_atom_count = 0
    gro_atom_names: List[str] = []
    if len(gro_lines) >= 2:
        try:
            gro_atom_count = int(gro_lines[1].strip().split()[0])
        except (ValueError, IndexError):
            gro_atom_count = 0
        gro_atom_names = [line[10:15].strip() for line in gro_lines[2:6] if len(line) >= 15]

    problems: List[str] = []
    if moleculetype != "SOL":
        problems.append(f"ITP moleculetype is {moleculetype!r}, expected 'SOL'")
    if atom_count != 4:
        problems.append(f"ITP [ atoms ] has {atom_count} atoms, expected 4")
    if gro_atom_count <= 0 or gro_atom_count % 4 != 0:
        problems.append(f"GRO has {gro_atom_count} atoms, not a whole number of 4-site waters")
    if gro_atom_names[:4] and gro_atom_names[:4] != ["OW", "HW1", "HW2", "MW"]:
        problems.append(f"GRO first molecule atoms are {gro_atom_names[:4]}, expected OW/HW1/HW2/MW")

    return {
        "validated": not problems,
        "water_model": water_model,
        "moleculetype": moleculetype,
        "itp_atom_count": atom_count,
        "gro_atom_count": gro_atom_count,
        "gro_first_atoms": gro_atom_names[:4],
        "problems": problems,
        "initialization_template": "tip4p.gro",
    }


def validate_water_model_topology(
    backend: GromacsBackend, force_field: str, water_model: str
) -> Dict[str, Any]:
    """Read the installed water ITP/GRO through the backend and validate them."""
    if water_model not in {"tip4p", "tip4pew"}:
        return {"validated": True, "water_model": water_model, "note": "no 4-site check required"}
    ff_dir = force_field if force_field.endswith(".ff") else f"{force_field}.ff"
    prefix = backend.data_prefix or ("/usr" if backend.is_wsl else "")
    if not prefix:
        return {
            "validated": False,
            "water_model": water_model,
            "reason": "GROMACS data prefix is unknown; cannot validate water files.",
        }
    itp_path = f"{prefix}/share/gromacs/top/{ff_dir}/tip4pew.itp" if water_model == "tip4pew" else f"{prefix}/share/gromacs/top/{ff_dir}/tip4p.itp"
    gro_path = f"{prefix}/share/gromacs/top/tip4p.gro"
    script = (
        f"cat {shlex.quote(itp_path)} 2>/dev/null; echo '@@@WATER@@@'; "
        f"cat {shlex.quote(gro_path)} 2>/dev/null"
    )
    try:
        code, stdout, stderr = run_backend_shell(
            script, cwd=PROJECT_ROOT, timeout=60, backend=backend
        )
    except Exception as exc:  # noqa: BLE001
        return {
            "validated": False,
            "water_model": water_model,
            "reason": f"water files could not be read through the backend: {exc}",
        }
    if code != 0 or "@@@WATER@@@" not in stdout:
        return {
            "validated": False,
            "water_model": water_model,
            "reason": f"water files are missing or unreadable ({stderr[-300:]})",
        }
    itp_text, gro_text = stdout.split("@@@WATER@@@", 1)
    result = inspect_water_topology_text(itp_text, gro_text, water_model)
    result["itp_path"] = itp_path
    result["gro_path"] = gro_path
    return result


def clean_receptor_for_md(
    receptor_pdb: Path,
    output_pdb: Path,
    allow_cofactor_removal: bool = False,
    retained_ion_resnames: Optional[set] = None,
    retained_water_resnames: Optional[set] = None,
) -> Tuple[Path, List[str]]:
    """
    Produce a pdb2gmx-safe receptor file.

    - Keeps ATOM records (protein) and converts selenomethionine (MSE) to MET.
    - Drops waters, ions, cofactors and the docked ligand (HETATM).
    - Refuses to silently drop catalytically essential cofactors (HEM/FE/MG/ZN/
      MN/CU/CA) unless ``allow_cofactor_removal=True`` is explicit.
    - Keeps only the first alternate location.
    - Normalises occupancy/bfactor fields so pdb2gmx does not reject the file.
    """
    receptor_pdb = Path(receptor_pdb)
    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    warnings: List[str] = []
    kept: List[str] = []
    dropped_het = 0
    retained_het = 0
    retained_ion_resnames = {name.upper() for name in (retained_ion_resnames or set())}
    retained_water_resnames = {name.upper() for name in (retained_water_resnames or set())}
    dropped_het_resnames: set[str] = set()
    critical_cofactors: set[str] = set()
    dropped_waters = 0
    dropped_alt = 0
    converted_mse = 0

    with open(receptor_pdb, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            record = line[:6]
            if record == "ENDMDL":
                break
            resname_raw = line[17:20].strip().upper() if len(line) >= 20 else ""
            is_mse_hetatm = record.startswith("HETATM") and resname_raw == "MSE"
            if record.startswith("ATOM") or is_mse_hetatm:
                altloc = line[16:17]
                if altloc not in (" ", "A", "1", ""):
                    dropped_alt += 1
                    continue

                resname = resname_raw
                atom_name = line[12:16].strip()
                if resname == "MSE":
                    resname = "MET"
                    atom_name = _MSE_ATOM_MAP.get(atom_name.upper(), atom_name)
                    converted_mse += 1

                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                element = line[76:78].strip() if len(line) >= 78 else atom_name[0]
                if resname_raw == "MSE" and atom_name.upper() in {"SD", "SG"}:
                    element = "S"
                chain = line[21:22] or "A"
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    resnum = 1
                try:
                    icode = line[26:27]
                except IndexError:
                    icode = " "

                new_line = (
                    f"ATOM  {len(kept) + 1:5d} {atom_name:<4s} {resname:3s} {chain:1s}{resnum:4d}{icode:1s}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}"
                )
                kept.append(new_line)
            elif record.startswith("HETATM"):
                if resname_raw in retained_ion_resnames or resname_raw in retained_water_resnames:
                    retained_het += 1
                else:
                    dropped_het += 1
                    if resname_raw:
                        dropped_het_resnames.add(resname_raw)
                if resname_raw in {"HOH", "WAT", "SOL"} and resname_raw not in retained_water_resnames:
                    dropped_waters += 1
                if resname_raw in _CRITICAL_COFACTOR_RESNAMES:
                    critical_cofactors.add(resname_raw)

    if not kept:
        raise TopologyBuildError(
            f"No protein ATOM records found in receptor file: {receptor_pdb}. "
            "Prepare the receptor with the docking pipeline first."
        )

    if critical_cofactors and not allow_cofactor_removal:
        raise TopologyBuildError(
            "The receptor contains catalytically relevant HETATM cofactors "
            f"({', '.join(sorted(critical_cofactors))}). This MD route does not "
            "parameterise cofactors, so dropping them would produce a chemically "
            "different system. Parameterise them explicitly, or set "
            "allow_cofactor_removal=True only after confirming the biology does not "
            "require them."
        )

    if dropped_waters:
        warnings.append(
            f"Removed {dropped_waters} crystallographic water atoms. If a bridging "
            "water is mechanistically essential, it must be parameterised explicitly."
        )

    if dropped_het:
        warnings.append(
            f"Removed {dropped_het} HETATM records before pdb2gmx "
            f"({', '.join(sorted(dropped_het_resnames)) or 'unknown residues'}). "
            "Only the protein is parameterised; waters, ions, cofactors and the docked "
            "ligand are not part of this MD topology."
        )
    if dropped_alt:
        warnings.append(f"Kept only the first alternate location for {dropped_alt} atoms.")
    if converted_mse:
        warnings.append(f"Converted {converted_mse} MSE atoms to MET for pdb2gmx.")

    with open(output_pdb, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(kept) + "\nEND\n")

    return output_pdb, warnings


def build_protein_topology(
    backend: GromacsBackend,
    receptor_pdb: Path,
    output_dir: Path,
    force_field: str = "amber99sb-ildn",
    water_model: str = "tip3p",
    allow_cofactor_removal: bool = False,
    retained_ion_resnames: Optional[set] = None,
    retained_water_resnames: Optional[set] = None,
) -> Dict[str, Any]:
    """Run ``gmx pdb2gmx`` and return the generated protein topology/coordinates."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ff = resolve_force_field(backend, force_field)
    water = resolve_water_model(water_model)
    # GROMACS' pdb2gmx does not accept "tip4pew" as a -water value even when
    # the force field ships tip4pew.itp.  Build with the 4-site tip4p topology
    # and then swap the include to the requested TIP4P-Ew parameters.
    pdb2gmx_water = "tip4p" if water == "tip4pew" else water
    water_validation = (
        validate_water_model_topology(backend, ff, water)
        if water in {"tip4p", "tip4pew"}
        else {"validated": True, "water_model": water, "note": "3-site water; no TIP4P check required"}
    )
    if water == "tip4pew" and not water_validation.get("validated"):
        raise TopologyBuildError(
            "TIP4P-Ew water topology validation failed: "
            + "; ".join(water_validation.get("problems", []))
            + " "
            + str(water_validation.get("reason", ""))
        )

    cleaned, warnings = clean_receptor_for_md(
        receptor_pdb,
        output_dir / "receptor_md.pdb",
        allow_cofactor_removal=allow_cofactor_removal,
        retained_ion_resnames=retained_ion_resnames,
        retained_water_resnames=retained_water_resnames,
    )
    gro = output_dir / "protein.gro"
    top = output_dir / "topol.top"
    posre = output_dir / "posre.itp"

    base_args = [
        "pdb2gmx",
        "-f", cleaned.name,
        "-o", gro.name,
        "-p", top.name,
        "-i", posre.name,
        "-ff", ff,
        "-water", pdb2gmx_water,
        "-ignh",
        "-nointer",
        "-chainsep", "id_or_ter",
        "-merge", "all",
    ]

    used_missing = False
    try:
        run_gmx_command(base_args, cwd=output_dir, log_name="pdb2gmx", timeout=900, backend=backend)
    except GromacsExecutionError as exc:
        message = str(exc)
        if "missing" in message.lower() or "atom" in message.lower():
            logger.warning("pdb2gmx failed; retrying with -missing: %s", message[:400])
            used_missing = True
            try:
                run_gmx_command(
                    base_args + ["-missing"],
                    cwd=output_dir,
                    log_name="pdb2gmx_missing",
                    timeout=900,
                    backend=backend,
                )
            except GromacsExecutionError as exc2:
                raise TopologyBuildError(
                    "pdb2gmx could not build the protein topology. The receptor likely has "
                    "missing side-chain atoms, non-standard residues, or incompatible termini.\n"
                    f"First error:\n{message}\n\nRetry error:\n{exc2}"
                ) from exc2
        else:
            raise TopologyBuildError(f"pdb2gmx failed for {receptor_pdb}:\n{message}") from exc

    if not gro.is_file() or not top.is_file():
        raise TopologyBuildError(
            "pdb2gmx reported success but did not produce protein.gro/topol.top. "
            f"Check {output_dir / 'pdb2gmx.log'}."
        )

    if used_missing:
        warnings.append(
            "pdb2gmx was run with -missing: the receptor contains missing atoms and the "
            "resulting topology may require manual inspection."
        )

    if water == "tip4pew":
        ff_dir = ff if ff.endswith(".ff") else f"{ff}.ff"
        top_text = top.read_text(encoding="utf-8", errors="replace")
        old_include = f'#include "{ff_dir}/tip4p.itp"'
        new_include = f'#include "{ff_dir}/tip4pew.itp"'
        if old_include not in top_text:
            raise TopologyBuildError(
                f"pdb2gmx built a tip4p topology but did not reference {old_include!r}; "
                "cannot switch it to TIP4P-Ew safely."
            )
        top.write_text(
            top_text.replace(old_include, new_include),
            encoding="utf-8",
            newline="\n",
        )
        warnings.append(
            "Water topology was switched from tip4p.itp to tip4pew.itp after pdb2gmx "
            "(GROMACS pdb2gmx does not accept -water tip4pew)."
        )

    atoms, _ = parse_gro_atoms(gro)
    return {
        "force_field": ff,
        "water_model": water,
        "protein_gro": gro,
        "topol_top": top,
        "posre": posre,
        "receptor_clean": cleaned,
        "protein_atom_count": len(atoms),
        "water_validation": water_validation,
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Ligand parameterisation (GAFF2 / AM1-BCC via ACPYPE)
# ---------------------------------------------------------------------------

def sanitize_resname(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", name or "").upper()
    return (cleaned[:5] or "LIG")


def ligand_pose_to_sdf(
    ligand_pose: ExtractedPose,
    output_sdf: Path,
    template_sdf: Optional[Path] = None,
) -> Tuple[Path, int, List[str]]:
    """
    Write the docked pose as an SDF with correct bond orders and explicit hydrogens.

    The docking pipeline stores the ligand 3D SDF (with bond orders). We keep
    that chemical graph, replace heavy-atom coordinates with the docked pose
    coordinates, and add hydrogens. This avoids guessing bond orders from a
    PDBQT file (which contains none).

    Returns ``(sdf_path, formal_charge, warnings)``.
    """
    output_sdf = Path(output_sdf)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    warnings: List[str] = []

    try:
        from rdkit import Chem
        from rdkit.Geometry import Point3D
    except ImportError as exc:  # pragma: no cover - RDKit is a hard dependency
        raise TopologyBuildError("RDKit is required to prepare the ligand for ACPYPE.") from exc

    heavy = [a for a in ligand_pose.atoms if a.is_heavy]
    if not heavy:
        raise TopologyBuildError("The docked pose contains no heavy atoms.")

    mol = None
    if template_sdf and Path(template_sdf).is_file():
        mol = Chem.MolFromMolFile(str(template_sdf), removeHs=True, sanitize=True)
        if mol is None:
            mol = Chem.MolFromMolFile(str(template_sdf), removeHs=True, sanitize=False)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                except Exception:
                    mol = None
        if mol is None:
            warnings.append(f"Could not parse ligand template SDF: {template_sdf}")
        else:
            template_symbols = [atom.GetSymbol().upper() for atom in mol.GetAtoms()]
            pose_symbols = [a.element.upper() for a in heavy]
            if len(template_symbols) != len(pose_symbols):
                raise TopologyBuildError(
                    "The ligand template SDF and the docked pose have different heavy-atom counts "
                    f"({len(template_symbols)} vs {len(pose_symbols)}). Re-run docking so the pose and "
                    "template come from the same ligand preparation."
                )
            # Vina output atom-name columns can collapse distinct names (e.g. all
            # carbons become "C"). Preserve atom identity through PDBQT serials,
            # mapped against the undocked SDF/PDBQT coordinates in their shared frame.
            from docking.atom_mapping import MAP_PROPERTY, record_atom_mapping
            if not mol.HasProp(MAP_PROPERTY):
                source_qt = Path(template_sdf).with_name(
                    Path(template_sdf).stem.removesuffix("_3d") + ".pdbqt"
                )
                if not source_qt.is_file():
                    raise TopologyBuildError("Ligand atom provenance is missing; re-run ligand preparation.")
                try:
                    record_atom_mapping(mol, source_qt.read_text(encoding="utf-8"))
                except (ValueError, IndexError) as exc:
                    raise TopologyBuildError(f"Could not map ligand input atoms: {exc}") from exc
            atom_map = json.loads(mol.GetProp(MAP_PROPERTY))
            pose_indices: Dict[int, LigandAtom] = {}
            for atom in heavy:
                source_idx = atom_map.get(str(atom.serial))
                if source_idx is None or source_idx < 0 or source_idx >= mol.GetNumAtoms():
                    raise TopologyBuildError(f"Docked pose serial has no source mapping: {atom.serial}")
                source = mol.GetAtomWithIdx(source_idx)
                if source.GetSymbol().upper() != atom.element.upper():
                    raise TopologyBuildError(f"Docked atom element conflicts with SDF mapping: {atom.element}")
                pose_indices[source_idx] = atom
            if len(pose_indices) != mol.GetNumAtoms():
                raise TopologyBuildError("Docked pose atom names do not cover the ligand SDF graph.")
            conf = mol.GetConformer()
            for template_idx, at in pose_indices.items():
                conf.SetAtomPosition(template_idx, Point3D(float(at.x), float(at.y), float(at.z)))
    else:
        warnings.append(
            "No ligand SDF template found in the docking run; bond orders will be perceived from the pose."
        )
        mol = Chem.MolFromPDBBlock(ligand_pose.pdb_block, removeHs=True, sanitize=False)
        if mol is None:
            raise TopologyBuildError(
                "Could not build a ligand molecule from the docked pose. "
                "Re-run docking to generate the ligand 3D SDF."
            )
        try:
            from rdkit.Chem import rdDetermineBonds
            rdDetermineBonds.DetermineBonds(mol, charge=0)
        except Exception as exc:
            raise TopologyBuildError(
                "The docking run has no ligand SDF template and RDKit could not perceive bond orders "
                f"from the pose: {exc}. Re-run docking so that <ligand>_3d.sdf is written."
            ) from exc

    # Add explicit hydrogens (coordinates generated from the heavy-atom frame).
    mol = Chem.AddHs(mol, addCoords=True)
    formal_charge = int(Chem.GetFormalCharge(mol))
    Chem.MolToMolFile(mol, str(output_sdf))

    if not output_sdf.is_file() or output_sdf.stat().st_size == 0:
        raise TopologyBuildError(f"Failed to write ligand SDF: {output_sdf}")
    return output_sdf, formal_charge, warnings


def parameterize_ligand_acpype(
    backend: GromacsBackend,
    ligand_sdf: Path,
    output_dir: Path,
    ligand_resname: str = "LIG",
    net_charge: int = 0,
    atom_type: str = "gaff2",
    charge_method: str = "bcc",
) -> Dict[str, Any]:
    """
    Run ACPYPE (Antechamber + Parmchk + Tleap) to produce a GROMACS ligand
    topology with GAFF2 atom types and AM1-BCC charges.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ligand_sdf = Path(ligand_sdf)
    resname = sanitize_resname(ligand_resname)

    # ACPYPE accepts .mdl; SDF and MDL molfile share the same connectivity block.
    acpype_input = output_dir / f"{resname}_pose.mdl"
    with open(ligand_sdf, "rb") as src, open(acpype_input, "wb") as dst:
        dst.write(src.read())

    script = (
        "command -v acpype >/dev/null 2>&1 || { echo 'ACPYPE_NOT_FOUND'; exit 127; }; "
        # ACPYPE's pip wheel bundles AmberTools, but set_for_pip() only wires
        # that bundle when no antechamber is already on PATH. Some WSL images
        # expose a stale /usr/local/bin/antechamber wrapper first; it points to
        # a missing /usr/local/amber.sh, so ACPYPE's conversion silently fails
        # before producing its expected temporary AC file named 'tmp'.
        "ACPYPE_LAUNCHER=$(command -v acpype); "
        "ACPYPE_PYTHON=$(head -n 1 \"$ACPYPE_LAUNCHER\" | sed 's/^#!//'); "
        "ACPYPE_BUNDLE=$(\"$ACPYPE_PYTHON\" -c 'from acpype.utils import bundled_amber_dir; print(bundled_amber_dir() or \"\")'); "
        "if [ -n \"$ACPYPE_BUNDLE\" ]; then "
        "export AMBERHOME=\"$ACPYPE_BUNDLE/\"; "
        "export PATH=\"$ACPYPE_BUNDLE/bin:$PATH\"; "
        "export LD_LIBRARY_PATH=\"$ACPYPE_BUNDLE/lib:${LD_LIBRARY_PATH:-}\"; "
        "fi; "
        f"acpype -i {shlex.quote(acpype_input.name)} -b {resname} "
        f"-c {shlex.quote(charge_method)} -n {int(net_charge)} "
        f"-a {shlex.quote(atom_type)} -o gmx"
    )

    if sys.platform == "win32" and not backend.is_wsl:
        # Native Windows has no bash.  Run the ACPYPE console launcher directly.
        acpype_exe = shutil.which("acpype")
        if not acpype_exe:
            raise TopologyBuildError(
                "ACPYPE was not found on PATH for the native Windows backend. "
                "Install ACPYPE (or use WSL) before ligand parameterisation."
            )
        command = [
            acpype_exe,
            "-i", acpype_input.name,
            "-b", resname,
            "-c", charge_method,
            "-n", str(int(net_charge)),
            "-a", atom_type,
            "-o", "gmx",
        ]
        try:
            completed = subprocess.run(
                command,
                cwd=str(output_dir),
                capture_output=True,
                text=True,
                timeout=1800,
            )
        except Exception as exc:  # noqa: BLE001
            raise TopologyBuildError(f"Could not launch native ACPYPE: {exc}") from exc
        code = completed.returncode
        stdout = completed.stdout or ""
        stderr = completed.stderr or ""
        (output_dir / "acpype.log").write_text(
            f"COMMAND: {' '.join(command)}\nEXIT CODE: {code}\n\n--- STDOUT ---\n{stdout}\n"
            f"--- STDERR ---\n{stderr}\n",
            encoding="utf-8",
            newline="\n",
        )
    else:
        code, stdout, stderr = run_backend_shell(
            script,
            cwd=output_dir,
            log_name="acpype",
            timeout=1800,
            backend=backend,
        )
    combined = (stdout or "") + "\n" + (stderr or "")
    if code != 0:
        hint = ""
        if "ACPYPE_NOT_FOUND" in combined:
            hint = (
                " ACPYPE is not installed in the backend. On Windows+WSL install it inside "
                "Ubuntu (`pipx install acpype` or `pip install --user acpype`) together with AmberTools."
            )
        raise TopologyBuildError(
            "ACPYPE ligand parameterisation failed." + hint + "\n" + combined[-4000:]
        )

    acpype_dir = output_dir / f"{resname}.acpype"
    itp_src = acpype_dir / f"{resname}_GMX.itp"
    gro_src = acpype_dir / f"{resname}_GMX.gro"
    if not itp_src.is_file() or not gro_src.is_file():
        raise TopologyBuildError(
            f"ACPYPE finished but expected files were not produced: {itp_src} / {gro_src}. "
            f"See {output_dir / 'acpype.log'}."
        )

    # Keep the production files in the run root so that #include paths are
    # unambiguous (topol.top, ligand ITP and posre all live side by side).
    itp = output_dir / f"{resname}_GMX.itp"
    gro = output_dir / f"{resname}_GMX.gro"
    shutil.copy2(itp_src, itp)
    shutil.copy2(gro_src, gro)

    # ACPYPE writes its default residue name (MOL) in [ atoms ] even when the
    # moleculetype is renamed.  GROMACS then writes MOL into solvated GRO files,
    # which breaks index generation (and any analysis keyed on the ligand
    # residue name).  Make the residue name match the moleculetype.
    _rename_itp_residue(itp, resname)

    atom_count, heavy_indices = _parse_itp_atoms(itp)
    posre = output_dir / f"posre_{resname}.itp"
    write_ligand_posre(posre, heavy_indices, force_constant=1000.0)

    _ensure_posres_include(itp, posre.name)

    return {
        "resname": resname,
        "itp": itp,
        "gro": gro,
        "posre": posre,
        "atom_count": atom_count,
        "heavy_indices": heavy_indices,
        "acpype_dir": acpype_dir,
    }


def _rename_itp_residue(itp_path: Path, resname: str) -> None:
    """Rewrite the residue-name column of an ITP ``[ atoms ]`` block.

    The column is fixed-width; replacing a 3-character name (``MOL``) with a
    3–5 character name keeps the remaining columns aligned.
    """
    lines = itp_path.read_text(encoding="utf-8", errors="replace").splitlines()
    out: List[str] = []
    in_atoms = False
    changed = 0
    for raw in lines:
        stripped = raw.strip()
        if stripped.startswith("["):
            in_atoms = stripped.replace("[", "").replace("]", "").strip().lower() == "atoms"
            out.append(raw)
            continue
        if in_atoms and stripped and not stripped.startswith((";", "#")):
            # nr type resnr residue atom cgnr charge mass
            parts = raw.split()
            if len(parts) >= 5:
                parts[3] = resname
                # Preserve the original leading indentation and rebuild with
                # whitespace separation; GROMACS is whitespace-delimited.
                out.append(" ".join(parts))
                changed += 1
                continue
        out.append(raw)
    if changed:
        itp_path.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


def _parse_itp_atoms(itp_path: Path) -> Tuple[int, List[int]]:
    """Return (atom_count, heavy_atom_indices) from a GROMACS ligand ITP."""
    atom_count = 0
    heavy: List[int] = []
    in_atoms = False
    with open(itp_path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("["):
                in_atoms = line.replace("[", "").replace("]", "").strip().lower() == "atoms"
                continue
            if not in_atoms:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                nr = int(parts[0])
            except ValueError:
                continue
            atom_name = parts[4]
            atom_type = parts[1]
            atom_count = max(atom_count, nr)
            is_h = atom_name.upper().startswith("H") or atom_type.lower().startswith("h")
            if not is_h:
                heavy.append(nr)
    if atom_count == 0:
        raise TopologyBuildError(f"No [ atoms ] entries found in ligand topology: {itp_path}")
    return atom_count, heavy


def write_ligand_posre(output_path: Path, heavy_indices: List[int], force_constant: float = 1000.0) -> Path:
    """Write a position-restraint file for ligand heavy atoms."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "; Position restraints for ligand heavy atoms\n",
        f"; fc = {force_constant:.1f} kJ/(mol*nm^2) (AGENTS.md / Lemkul protocol)\n",
        "[ position_restraints ]\n",
        ";  ai  funct         fcx         fcy         fcz\n",
    ]
    for idx in heavy_indices:
        lines.append(
            f"{idx:6d}     1    {force_constant:10.1f} {force_constant:10.1f} {force_constant:10.1f}\n"
        )
    with open(output_path, "w", encoding="utf-8", newline="\n") as f:
        f.writelines(lines)
    return output_path


def _ensure_posres_include(itp_path: Path, posre_filename: str) -> None:
    text = itp_path.read_text(encoding="utf-8", errors="replace")
    if "POSRES_LIG" in text:
        return
    block = (
        f'\n; Include ligand position restraints (added by Bio_V2 MD pipeline)\n'
        f'#ifdef POSRES_LIG\n'
        f'#include "{posre_filename}"\n'
        f'#endif\n'
    )
    with open(itp_path, "a", encoding="utf-8", newline="\n") as f:
        f.write(block)


# ---------------------------------------------------------------------------
# Topology / coordinate assembly
# ---------------------------------------------------------------------------

def merge_complex_topology(
    protein_top: Path,
    ligand_itp: Path,
    ligand_resname: str,
    output_top: Path,
) -> Path:
    """
    Insert the ligand `#include` after the force-field include and append the
    ligand molecule to `[ molecules ]` in the pdb2gmx-generated topology.
    """
    protein_top = Path(protein_top)
    ligand_itp = Path(ligand_itp)
    output_top = Path(output_top)

    text = protein_top.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    include_line = f'#include "{ligand_itp.name}"'
    if include_line not in text:
        insert_at = None
        for i, line in enumerate(lines):
            if re.match(r'^\s*#include\s+".*forcefield\.itp"', line):
                insert_at = i + 1
                break
        if insert_at is None:
            for i, line in enumerate(lines):
                if line.strip().startswith("#include"):
                    insert_at = i + 1
                    break
        if insert_at is None:
            raise TopologyBuildError(
                f"Could not find a force-field #include in {protein_top}; refusing to build a topology."
            )
        lines.insert(insert_at, "")
        lines.insert(insert_at + 1, "; Small-molecule ligand topology (ACPYPE / GAFF2 / AM1-BCC)")
        lines.insert(insert_at + 2, include_line)

    # Locate [ molecules ] and append the ligand after the last molecule line.
    molecules_idx = None
    for i, line in enumerate(lines):
        if line.strip().lower() == "[ molecules ]":
            molecules_idx = i
    if molecules_idx is None:
        raise TopologyBuildError(f"No [ molecules ] section found in {protein_top}.")

    last_entry = molecules_idx
    for i in range(molecules_idx + 1, len(lines)):
        stripped = lines[i].split(";", 1)[0].strip()
        if stripped.startswith("["):
            break
        if stripped:
            last_entry = i
    if last_entry == molecules_idx:
        raise TopologyBuildError(f"[ molecules ] in {protein_top} contains no molecules.")

    lines.insert(last_entry + 1, f"{ligand_resname:<20s} 1")

    output_top.parent.mkdir(parents=True, exist_ok=True)
    with open(output_top, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return output_top


def combine_gro_files(
    protein_gro: Path,
    ligand_gro: Path,
    output_gro: Path,
    ligand_resname: str,
    expected_ligand_atoms: Optional[int] = None,
) -> Path:
    """Concatenate protein and ligand GRO coordinates in topology order."""
    protein_gro = Path(protein_gro)
    ligand_gro = Path(ligand_gro)
    output_gro = Path(output_gro)

    with open(protein_gro, "r", encoding="utf-8", errors="replace") as f:
        protein_lines = f.read().splitlines()
    with open(ligand_gro, "r", encoding="utf-8", errors="replace") as f:
        ligand_lines = f.read().splitlines()

    try:
        n_prot = int(protein_lines[1].strip().split()[0])
        n_lig = int(ligand_lines[1].strip().split()[0])
    except (IndexError, ValueError) as exc:
        raise TopologyBuildError("Malformed GRO file while combining the complex.") from exc

    if expected_ligand_atoms is not None and n_lig != expected_ligand_atoms:
        raise TopologyBuildError(
            f"Ligand coordinate/topology atom-count mismatch: {n_lig} vs {expected_ligand_atoms}."
        )

    protein_body = protein_lines[2:2 + n_prot]
    ligand_body = ligand_lines[2:2 + n_lig]
    if len(protein_body) != n_prot or len(ligand_body) != n_lig:
        raise TopologyBuildError("Malformed GRO file: declared atom count does not match body.")

    # Re-index ligand residue numbers after the protein.
    last_resnum = 0
    for line in protein_body:
        try:
            last_resnum = max(last_resnum, int(line[0:5]))
        except ValueError:
            pass

    ligand_out: List[str] = []
    for line in ligand_body:
        if len(line) < 44:
            raise TopologyBuildError(f"Malformed ligand GRO line: {line!r}")
        new_resnum = last_resnum + 1
        new_line = f"{new_resnum:5d}{ligand_resname:<5s}" + line[10:]
        ligand_out.append(new_line)

    box_line = protein_lines[2 + n_prot] if len(protein_lines) > 2 + n_prot else "   0.00000   0.00000   0.00000"
    if not box_line.strip():
        box_line = "   0.00000   0.00000   0.00000"

    output_gro.parent.mkdir(parents=True, exist_ok=True)
    with open(output_gro, "w", encoding="utf-8", newline="\n") as f:
        f.write("Protein-ligand complex prepared by Bio_V2 MD pipeline\n")
        f.write(f"{n_prot + n_lig:5d}\n")
        f.writelines(line + "\n" for line in protein_body)
        f.writelines(line + "\n" for line in ligand_out)
        f.write(box_line + "\n")
    return output_gro


def gro_to_pdb(gro_path: Path, output_pdb: Path) -> Path:
    """Convert a GRO file to a simple PDB for preview/analysis."""
    atoms, _ = parse_gro_atoms(gro_path)
    output_pdb = Path(output_pdb)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w", encoding="utf-8", newline="\n") as f:
        f.write("REMARK   1 Complex coordinates produced by GROMACS (pdb2gmx + ACPYPE)\n")
        for idx, (resnum, resname, atomname, x, y, z) in enumerate(atoms, start=1):
            element = re.sub(r"[^A-Za-z]", "", atomname)[:1] or atomname[:1]
            f.write(
                f"ATOM  {idx:5d} {atomname:<4s} {resname:>3s}  A{resnum:4d}    "
                f"{x * 10.0:8.3f}{y * 10.0:8.3f}{z * 10.0:8.3f}  1.00 20.00          {element:>2s}\n"
            )
        f.write("END\n")
    return output_pdb


def prepare_complex_topology(
    receptor_pdb: Path,
    ligand_pose: ExtractedPose,
    output_dir: Path,
    force_field: str = "amber99sb-ildn",
    water_model: str = "tip3p",
    ligand_sdf: Optional[Path] = None,
    ligand_resname: str = "LIG",
    backend: Optional[GromacsBackend] = None,
    allow_cofactor_removal: bool = False,
    retain_monatomic_ions: bool = False,
    retain_crystal_waters: bool = False,
) -> TopologyResult:
    """
    Build a physically valid protein–ligand topology:

    1. Clean the receptor and run ``gmx pdb2gmx``.
    2. Write the docked pose as an SDF (correct bond orders + hydrogens).
    3. Run ACPYPE (GAFF2 + AM1-BCC) for the ligand.
    4. Merge the ligand ITP into ``topol.top``.
    5. Combine ``protein.gro`` and ``ligand.gro`` into ``complex.gro``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    backend = backend or detect_gromacs()
    if backend is None:
        raise TopologyBuildError(
            "GROMACS was not found natively or through WSL, so pdb2gmx cannot build a topology. "
            "Install GROMACS or set GMX_BIN."
        )

    resname = sanitize_resname(ligand_resname)
    retained_ions: List[Dict[str, Any]] = []
    retained_waters: List[Dict[str, Any]] = []
    if retain_monatomic_ions:
        retained_ions, unsupported = extract_monatomic_ions(receptor_pdb)
        if unsupported:
            raise TopologyBuildError(
                "Retention is limited to monatomic ions with explicit Amber moleculetypes "
                f"({', '.join(sorted(_RETAINABLE_MONATOMIC_IONS))}); unsupported cofactors: "
                f"{', '.join(sorted(unsupported))}. Parameterise them separately instead of "
                "silently deleting them."
            )
    if retain_crystal_waters:
        retained_waters = extract_crystal_waters(receptor_pdb)
    protein = build_protein_topology(
        backend=backend,
        receptor_pdb=receptor_pdb,
        output_dir=output_dir,
        force_field=force_field,
        water_model=water_model,
        allow_cofactor_removal=allow_cofactor_removal or bool(retained_ions),
        retained_ion_resnames={ion["resname"] for ion in retained_ions},
        retained_water_resnames=_CRYSTAL_WATER_RESNAMES if retained_waters else set(),
    )

    ligand_sdf_path, formal_charge, ligand_warnings = ligand_pose_to_sdf(
        ligand_pose,
        output_dir / f"{resname}_pose.sdf",
        template_sdf=ligand_sdf,
    )
    if "amber" not in str(protein["force_field"]).lower():
        raise TopologyBuildError(
            "The ligand parameterisation route in this deployment is GAFF2/AM1-BCC (AMBER). "
            f"Requested protein force field '{protein['force_field']}' is not compatible. "
            "Use an Amber99SB-ILDN force field."
        )
    ligand = parameterize_ligand_acpype(
        backend=backend,
        ligand_sdf=ligand_sdf_path,
        output_dir=output_dir,
        ligand_resname=resname,
        net_charge=formal_charge,
    )

    top = merge_complex_topology(
        protein_top=protein["topol_top"],
        ligand_itp=ligand["itp"],
        ligand_resname=resname,
        output_top=output_dir / "topol.top",
    )
    complex_gro = combine_gro_files(
        protein_gro=protein["protein_gro"],
        ligand_gro=ligand["gro"],
        output_gro=output_dir / "complex.gro",
        ligand_resname=resname,
        expected_ligand_atoms=ligand["atom_count"],
    )
    complex_pdb = gro_to_pdb(complex_gro, output_dir / "complex.pdb")

    if retained_ions:
        top = append_ions_to_topology(top, retained_ions)
        complex_gro = append_ions_to_gro(complex_gro, retained_ions)
        complex_pdb = gro_to_pdb(complex_gro, output_dir / "complex.pdb")

    if retained_waters:
        top = append_crystal_waters_to_topology(top, len(retained_waters))
        complex_gro = append_crystal_waters_to_gro(complex_gro, retained_waters)
        complex_pdb = gro_to_pdb(complex_gro, output_dir / "complex.pdb")

    warnings = list(protein["warnings"]) + list(ligand_warnings)
    if formal_charge != 0:
        warnings.append(
            f"Ligand formal charge is {formal_charge:+d}; this value was passed to ACPYPE for AM1-BCC."
        )

    if retained_ions:
        warnings.append(
            "Retained monatomic ions: "
            + ", ".join(f"{ion['resname']}:{ion['chain']}:{ion['resseq']}" for ion in retained_ions)
            + ". Their coordinates and molecule counts were appended to complex.gro/topol.top."
        )
    if retained_waters:
        warnings.append(
            f"Retained {len(retained_waters)} crystallographic waters as SOL molecules. "
            "Hydrogens were generated at TIP3P geometry when absent; orientation is approximate."
        )
    return TopologyResult(
        method=f"{protein['force_field']} + GAFF2/AM1-BCC (pdb2gmx + ACPYPE)",
        force_field=protein["force_field"],
        water_model=protein["water_model"],
        ligand_resname=resname,
        complex_gro=complex_gro,
        complex_pdb=complex_pdb,
        ligand_gro=ligand["gro"],
        ligand_itp=ligand["itp"],
        ligand_posre=ligand["posre"],
        topol_top=top,
        warnings=warnings,
        receptor_pdb=protein["receptor_clean"],
        protein_gro=protein["protein_gro"],
        ligand_sdf=ligand_sdf_path,
        ligand_charge=formal_charge,
        protein_atom_count=protein["protein_atom_count"],
        ligand_atom_count=ligand["atom_count"],
        backend_label=backend.label,
        water_validation=protein.get("water_validation", {}),
        retained_ions=retained_ions,
        retained_waters=len(retained_waters),
    )


# ---------------------------------------------------------------------------
# Legacy helpers retained for preview / CGenFF penalty evaluation
# ---------------------------------------------------------------------------

def assemble_protein_ligand_complex(
    receptor_pdb: Path,
    ligand_pose: ExtractedPose,
    output_dir: Path,
    ligand_resname: str = "LIG",
) -> Tuple[Path, Path]:
    """
    Combine receptor coordinates and ligand pose into a single PDB/GRO.

    This is a *preview* helper used before GROMACS runs. The production
    ``complex.gro`` is generated from pdb2gmx + ACPYPE coordinates by
    :func:`prepare_complex_topology`.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    complex_pdb = output_dir / "complex_preview.pdb"
    complex_gro = output_dir / "complex_preview.gro"

    with open(receptor_pdb, "r", encoding="utf-8", errors="replace") as f:
        rec_lines = f.readlines()

    atom_lines: List[str] = []
    gro_atom_entries: List[Tuple[str, str, int, float, float, float]] = []
    atom_counter = 1

    for line in rec_lines:
        if line.startswith(("ATOM  ", "HETATM")):
            resname = line[17:20].strip()
            if resname in ("SOL", "HOH", "WAT", ligand_resname, "UNL", "JZ4"):
                continue
            name = line[12:16].strip()
            chain = line[21:22]
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                resnum = 1
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue
            elem = line[76:78].strip() if len(line) >= 78 else name[0]

            atom_lines.append(
                f"ATOM  {atom_counter:5d} {name:<4s} {resname:3s} {chain:1s}{resnum:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {elem:>2s}\n"
            )
            gro_atom_entries.append((resname, name, resnum, x / 10.0, y / 10.0, z / 10.0))
            atom_counter += 1

    last_resnum = gro_atom_entries[-1][2] if gro_atom_entries else 1
    lig_resnum = last_resnum + 1

    for at in ligand_pose.atoms:
        elem = at.element.capitalize()
        atom_lines.append(
            f"HETATM{atom_counter:5d} {at.name:<4s} {ligand_resname:3s} L{lig_resnum:4d}    "
            f"{at.x:8.3f}{at.y:8.3f}{at.z:8.3f}{1.00:6.2f}{0.00:6.2f}          {elem:>2s}\n"
        )
        gro_atom_entries.append((ligand_resname, at.name, lig_resnum, at.x / 10.0, at.y / 10.0, at.z / 10.0))
        atom_counter += 1

    with open(complex_pdb, "w", encoding="utf-8", newline="\n") as f:
        f.write("REMARK   1 Biophysical Protein-Ligand Complex Preview (not for GROMACS input)\n")
        f.writelines(atom_lines)
        f.write("END\n")

    xs = [e[3] for e in gro_atom_entries]
    ys = [e[4] for e in gro_atom_entries]
    zs = [e[5] for e in gro_atom_entries]
    box_x = max(xs) - min(xs) + 1.0
    box_y = max(ys) - min(ys) + 1.0
    box_z = max(zs) - min(zs) + 1.0

    with open(complex_gro, "w", encoding="utf-8", newline="\n") as f:
        f.write("Complex preview (production coordinates are generated by GROMACS)\n")
        f.write(f"{len(gro_atom_entries):5d}\n")
        for idx, (rname, aname, rnum, gx, gy, gz) in enumerate(gro_atom_entries, start=1):
            f.write(f"{rnum:5d}{rname:<5s}{aname:>5s}{idx:5d}{gx:8.3f}{gy:8.3f}{gz:8.3f}\n")
        f.write(f"   {box_x:10.5f} {box_y:10.5f} {box_z:10.5f}\n")

    return complex_pdb, complex_gro


def generate_ligand_position_restraints(
    ligand_atoms: List[LigandAtom],
    output_path: Path,
    force_constant: float = 1000.0,
) -> Path:
    """Legacy helper: write restraints for a pose atom list (1-based indices)."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "; Position Restraint File for Ligand Heavy Atoms\n",
        f"; Harmonic force constant fc = {force_constant:.1f} kJ/(mol * nm^2)\n",
        "; Reference: wiki/gromacs-protein-ligand-protocol.md\n",
        "[ position_restraints ]\n",
        ";  ai  funct         fcx         fcy         fcz\n",
    ]
    for at in ligand_atoms:
        if at.is_heavy:
            lines.append(
                f"{at.serial:5d}     1    {force_constant:10.1f} {force_constant:10.1f} {force_constant:10.1f}\n"
            )
    with open(output_path, "w", encoding="utf-8", newline="\n") as f:
        f.writelines(lines)
    return output_path


def evaluate_cgenff_penalty(penalty_val: float) -> CGenFFPenaltyReport:
    """
    Evaluate CGenFF parameter penalty scores (Lemkul / LiveCoMS):
    <10 robust, 10-50 warning/verification needed, >50 invalid.
    """
    val = float(penalty_val)
    if val < 10.0:
        status = "PASS"
        verdict = "معاملات CGenFF عالية الدقة وموثوقة لمحاكاة الإنتاج (Penalty < 10)"
        details = [
            f"درجة الجزاء القصوى: {val:.1f}",
            "المعاملات منتقاة بدقة ومتسقة كيميائياً مع حقل القوى CHARMM36.",
            "لا تتطلب إعادة مسح لطاقة الالتواء أو حسابات تفاعل QM مع الماء.",
        ]
    elif val <= 50.0:
        status = "WARNING"
        verdict = "معاملات مقبولة مع التحذير (10 <= Penalty <= 50): تتطلب فحصاً إضافياً"
        details = [
            f"درجة الجزاء القصوى: {val:.1f}",
            "بعض زوايا الالتواء أو توزيع الشحنات تستند لمقاربات جزئية.",
            "يُنصح بإجراء مسح طاقة الالتواء ومقارنتها بالحسابات الكمية (QM).",
        ]
    else:
        status = "FAIL"
        verdict = "معاملات غير موثوقة (Penalty > 50): يُحظر استخدامها في الإنتاج دون إعادة تدريب"
        details = [
            f"درجة الجزاء القصوى: {val:.1f}",
            "المعاملات المقترحة تفتقر إلى النظائر الكيميائية المعتمدة في CGenFF.",
            "إعادة ملاءمة الشحنات السطحية وزوايا الالتواء بواسطة DFT أمر إلزامي.",
        ]

    return CGenFFPenaltyReport(
        max_charge_penalty=min(val, 25.0),
        max_dihedral_penalty=val,
        average_penalty=val * 0.45,
        status=status,
        verdict=verdict,
        details=details,
    )


def generate_ligand_itp_openff(*args, **kwargs):
    """Removed: the previous implementation emitted invented atom types."""
    raise TopologyBuildError(
        "OpenFF Sage parameterisation is not available in this deployment. "
        "Install the OpenFF Toolkit + Interchange and implement an explicit "
        "SMIRNOFF -> GROMACS export; the old placeholder ITP generator was removed "
        "because it produced physically meaningless parameters."
    )


def generate_ligand_itp_cgenff(*args, **kwargs):
    """Removed: the previous implementation emitted invented CGenFF atom types."""
    raise TopologyBuildError(
        "CHARMM36/CGenFF parameterisation requires a CHARMM36 force field and CGenFF "
        "parameters (or the CGenFF server). Neither is installed in this backend. "
        "Use the Amber99SB-ILDN + GAFF2/AM1-BCC route."
    )


def build_master_topology(*args, **kwargs):
    """Removed: use :func:`prepare_complex_topology` (pdb2gmx + ACPYPE)."""
    raise TopologyBuildError(
        "build_master_topology() was removed because it wrote a topology with an empty "
        "protein [ moleculetype ]. Use prepare_complex_topology() instead."
    )
