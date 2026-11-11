#!/usr/bin/env python3
"""convert_amber_to_gmx.py — deterministic ParmEd fallback for ACPYPE output.

Why this file exists
--------------------
ACPYPE occasionally aborts (``free(): invalid pointer``) *after* writing the
Amber files (``*_AC.prmtop`` / ``*_AC.inpcrd``) but *before* writing the
GROMACS files. This converter turns those Amber files into a GROMACS
``.itp``/``.gro`` pair without rerunning the expensive AM1-BCC calculation.

Scientific guard (PLAN.md Phase 2, wiki/force-fields.md)
--------------------------------------------------------
ParmEd copies charges verbatim, while ACPYPE calls ``balanceCharges()``
(acpype/topol.py:1564) to remove the rounding residual left by sqm. Without
that step the ligand carries e.g. +0.002 e, which makes the whole system
non-neutral and triggers GROMACS' PME background-charge warning. We therefore
apply the *same* algorithm as ACPYPE: distribute the residual equally over the
atoms that hold the extreme charge (max |q|).

Notes
-----
* The ``[ defaults ]`` section is still stripped: amber99sb-ildn.ff already
  declares ``comb-rule 2, fudgeLJ 0.5, fudgeQQ 0.8333`` (verified by
  scripts/verify_stack.py), so a duplicate block would be rejected.
* No silent success: every failure path exits non-zero with a clear message.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

import parmed as pmd

CHARGE_TOL = 1e-4


def itp_charge_sum(itp_path: Path) -> Tuple[float, int]:
    """Return (sum of [ atoms ] charges, atom count) from a GROMACS ITP."""
    total = 0.0
    count = 0
    in_atoms = False
    for raw in Path(itp_path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            in_atoms = line.replace("[", "").replace("]", "").strip().lower() == "atoms"
            continue
        if not in_atoms:
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            q = float(parts[6])
        except ValueError:
            continue
        total += q
        count += 1
    return total, count


def balance_charges(charges: List[float], target: float) -> Tuple[List[float], float]:
    """Mirror acpype.topol.Topology.balanceCharges().

    residual = total - target is removed from every atom whose |charge| equals
    the maximum |charge| (ACPYPE behaviour). Returns (new_charges, residual).
    """
    total = sum(charges)
    residual = total - target
    if abs(residual) < 1e-12:
        return list(charges), residual

    # Match ACPYPE exactly: choose the positive maximum on an absolute-value
    # tie; otherwise choose whichever signed extreme has larger magnitude.
    max_charge = max(charges)
    min_charge = min(charges)
    lim = max_charge if abs(max_charge) >= abs(min_charge) else min_charge
    n_lim = sum(1 for q in charges if q == lim)
    fixed = [lim - residual / n_lim if q == lim else q for q in charges]
    return fixed, residual


def strip_duplicate_sections(itp_path: Path) -> None:
    """Remove [ defaults ] / [ system ] / [ molecules ] and rename moleculetype."""
    lines = itp_path.read_text(encoding="utf-8").splitlines()
    clean: List[str] = []
    skip_section = False
    in_moltype = False
    resname = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("["):
            section = stripped.strip("[] \t").lower()
            if section in ("defaults", "system", "molecules"):
                skip_section = True
                in_moltype = False
                continue
            skip_section = False
            in_moltype = section == "moleculetype"
            clean.append(line)
            continue
        if skip_section:
            continue
        if in_moltype and stripped and not stripped.startswith(";"):
            parts = stripped.split()
            if len(parts) >= 2:
                resname = parts[0]
                clean.append(f" {parts[0]:<16s} {parts[1]}")
                in_moltype = False
                continue
        clean.append(line)
    itp_path.write_text("\n".join(clean) + "\n", encoding="utf-8", newline="\n")


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("resname", help="Ligand residue name, e.g. ERLOT")
    parser.add_argument(
        "--net-charge",
        type=int,
        default=None,
        help="Requested net charge (e). Default: round(sum of Amber charges), matching ACPYPE.",
    )
    parser.add_argument("--dir", default=None, help="ACPYPE output directory (default: <resname>.acpype)")
    args = parser.parse_args(argv)

    resname = args.resname
    base = Path(args.dir) if args.dir else Path(f"{resname}.acpype")
    prmtop = base / f"{resname}_AC.prmtop"
    inpcrd = base / f"{resname}_AC.inpcrd"
    itp_out = base / f"{resname}_GMX.itp"
    gro_out = base / f"{resname}_GMX.gro"

    if not (prmtop.is_file() and inpcrd.is_file()):
        print(f"ERROR: missing {prmtop} or {inpcrd}", file=sys.stderr)
        return 1

    mol = pmd.load_file(str(prmtop), str(inpcrd))
    for residue in mol.residues:
        residue.name = resname

    charges = [atom.charge for atom in mol.atoms]
    if args.net_charge is None:
        # ACPYPE: self.totalCharge = int(round(totalCharge / qConv)).
        target = float(int(round(sum(charges))))
        target_source = "auto-rounded"
    else:
        target = float(args.net_charge)
        target_source = "explicit"
    balanced, residual = balance_charges(charges, target)
    if abs(residual) > 1e-12:
        for atom, charge in zip(mol.atoms, balanced):
            atom.charge = charge

    mol.save(str(gro_out), overwrite=True)
    mol.save(str(itp_out), format="gromacs", overwrite=True)
    strip_duplicate_sections(itp_out)

    total, count = itp_charge_sum(itp_out)
    print(
        f"ParmEd fallback: atoms={count} residual_before={residual:+.8f} e "
        f"sum_after={total:+.8f} e target={target:+.1f} e ({target_source})"
    )
    if not itp_out.is_file() or not gro_out.is_file() or itp_out.stat().st_size == 0:
        print("ERROR: converter did not produce non-empty outputs", file=sys.stderr)
        return 1
    if abs(total - target) > CHARGE_TOL:
        print(
            f"ERROR: charge balance failed: sum={total:+.6f} e, "
            f"target={target:+.1f} e",
            file=sys.stderr,
        )
        return 1

    print(f"SUCCESS: Generated {itp_out} and {gro_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
