"""
GROMACS Index Generation from the Actual Solvated System
========================================================

Scientific Traceability & Wiki Documentation:
---------------------------------------------
- `wiki/gromacs-protein-ligand-protocol.md`: index groups for thermostat
  coupling (`Protein_Ligand`, `Water_and_ions`) and PBC post-processing
  (`Protein`, `Backbone`).
- `AGENTS.md`: merge protein and ligand into a single thermostat group.

Why this module exists
----------------------
Earlier versions wrote ``index.ndx`` from estimated atom counts (e.g. "4500
waters, 18 Na+, 12 Cl-"). That is scientifically invalid: the index must be
derived from the coordinates that GROMACS actually produced after
``solvate``/``genion``. This module parses a ``.gro`` file and builds the
groups deterministically.

The module is dependency-free and can be copied verbatim into an exported
GROMACS bundle and run as ``python make_index.py solv_ions.gro index.ndx LIG``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

WATER_RESNAMES = {"SOL", "HOH", "WAT", "TIP", "TIP3", "TIP4", "SPC"}
ION_RESNAMES = {"NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "BR", "IOD", "LI", "RB", "CS"}


def parse_gro_atoms(gro_path: Path) -> Tuple[List[Tuple[int, str, str, float, float, float]], str]:
    """
    Parse a GROMACS ``.gro`` file.

    Returns ``(atoms, box_line)`` where each atom tuple is
    ``(resnum, resname, atomname, x_nm, y_nm, z_nm)``.
    """
    gro_path = Path(gro_path)
    with open(gro_path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    if len(lines) < 3:
        raise ValueError(f"Malformed GRO file (too short): {gro_path}")

    try:
        declared = int(lines[1].strip().split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Malformed GRO atom-count line in {gro_path}") from exc

    body = lines[2:2 + declared]
    if len(body) != declared:
        raise ValueError(
            f"GRO file declares {declared} atoms but contains {len(body)} atom lines: {gro_path}"
        )

    atoms: List[Tuple[int, str, str, float, float, float]] = []
    for line in body:
        if len(line) < 44:
            raise ValueError(f"Malformed GRO atom line: {line!r}")
        try:
            resnum = int(line[0:5])
        except ValueError:
            resnum = 1
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        try:
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError as exc:
            raise ValueError(f"Malformed GRO coordinates: {line!r}") from exc
        atoms.append((resnum, resname, atomname, x, y, z))

    box_line = lines[2 + declared] if len(lines) > 2 + declared else ""
    return atoms, box_line


def _format_numbers(numbers: Iterable[int], cols: int = 15) -> List[str]:
    out: List[str] = []
    chunk: List[str] = []
    for num in numbers:
        chunk.append(f"{num:6d}")
        if len(chunk) >= cols:
            out.append(" ".join(chunk))
            chunk = []
    if chunk:
        out.append(" ".join(chunk))
    return out


def _is_hydrogen(atom_name: str) -> bool:
    # GROMACS atom names place the element first: H, H1, HA, HG1, 1H ...
    name = atom_name.strip()
    if not name:
        return False
    return name[0].upper() == "H" or (name[0].isdigit() and len(name) > 1 and name[1].upper() == "H")


def build_index_groups(
    gro_path: Path,
    ligand_resname: str = "LIG",
) -> Dict[str, List[int]]:
    """Build the required GROMACS index groups from the actual system."""
    atoms, _ = parse_gro_atoms(gro_path)
    lig = ligand_resname.strip().upper()

    protein: List[int] = []
    ligand: List[int] = []
    water_ions: List[int] = []
    water: List[int] = []
    ions: List[int] = []
    ligand_heavy: List[int] = []
    backbone: List[int] = []
    c_alpha: List[int] = []

    for idx, (_resnum, resname, atomname, _x, _y, _z) in enumerate(atoms, start=1):
        rn = resname.strip().upper()
        an = atomname.strip().upper()
        if rn == lig:
            ligand.append(idx)
            if not _is_hydrogen(atomname):
                ligand_heavy.append(idx)
        elif rn in WATER_RESNAMES:
            water.append(idx)
            water_ions.append(idx)
        elif rn in ION_RESNAMES:
            ions.append(idx)
            water_ions.append(idx)
        else:
            protein.append(idx)
            if an in {"N", "CA", "C", "O", "OXT", "OT1", "OT2"}:
                backbone.append(idx)
            if an == "CA":
                c_alpha.append(idx)

    if not protein:
        raise ValueError(f"No protein atoms detected in {gro_path}")
    if not ligand:
        raise ValueError(
            f"No ligand atoms with residue name '{ligand_resname}' detected in {gro_path}. "
            "Check that the ligand topology and coordinate residue names match."
        )

    groups: Dict[str, List[int]] = {
        "System": list(range(1, len(atoms) + 1)),
        "Protein": protein,
        "Ligand": ligand,
        "Protein_Ligand": protein + ligand,
        "Water_and_ions": water_ions,
        "Ligand_Heavy": ligand_heavy,
        "Backbone": backbone,
        "C-alpha": c_alpha,
    }
    if water:
        groups["Water"] = water
    if ions:
        groups["Ions"] = ions
    return groups


def write_index_file(groups: Dict[str, List[int]], output_ndx: Path) -> Path:
    output_ndx = Path(output_ndx)
    output_ndx.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    for name, members in groups.items():
        lines.append(f"[ {name} ]")
        lines.extend(_format_numbers(members))
    with open(output_ndx, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return output_ndx


def write_index_from_gro(
    gro_path: Path,
    output_ndx: Path,
    ligand_resname: str = "LIG",
) -> Dict[str, int]:
    """Build and write ``index.ndx``; returns group sizes."""
    groups = build_index_groups(gro_path, ligand_resname=ligand_resname)
    write_index_file(groups, output_ndx)
    return {name: len(members) for name, members in groups.items()}


def main(argv: Sequence[str]) -> int:
    if len(argv) < 3:
        print("Usage: python make_index.py <solv_ions.gro> <index.ndx> [ligand_resname]")
        return 2
    gro = Path(argv[1])
    ndx = Path(argv[2])
    lig = argv[3] if len(argv) > 3 else "LIG"
    sizes = write_index_from_gro(gro, ndx, ligand_resname=lig)
    print(f"Wrote {ndx} from {gro}")
    for name, size in sizes.items():
        print(f"  [ {name} ] = {size} atoms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
