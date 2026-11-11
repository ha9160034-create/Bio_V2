#!/opt/conda/envs/amber/bin/python
"""acpype_wrapper.py - deterministic ACPYPE launcher with a verified ParmEd fallback.

Why this wrapper exists
-----------------------
ACPYPE can abort with ``free(): invalid pointer`` *after* writing the Amber
files (``<root>_AC.prmtop`` / ``<root>_AC.inpcrd``) but *before* writing the
GROMACS pair. Re-running AM1-BCC is expensive, so we convert the already
computed Amber files with ParmEd instead.

Design references
-----------------
* PLAN.md Phase 2: four-state wrapper, charge neutrality, no silent success.
* wiki/open-force-field-and-ligand-parameterization.md: GAFF2 + AM1-BCC.
* wiki/force-fields.md: amber99sb-ildn / GAFF2 compatibility.
* wiki/gromacs-protein-ligand-protocol.md: ligand ITP/GRO contract.

Four-state policy (never overwrite a valid ACPYPE result)
---------------------------------------------------------
    ACPYPE rc | Amber files | GMX files             | action
    ----------+-------------+-----------------------+-------------------------
        0     | any         | valid + neutral       | accept as-is (no rewrite)
       !=0     | any         | valid + neutral       | accept as-is (rc is logged)
      any     | complete    | missing/bad/unbalanced| ParmEd convert + validate
      any     | missing     | missing/bad/unbalanced| fail (3 / 4 / 5)

The decision is made from the files, exactly as PLAN.md 2.1 specifies: a
non-zero ACPYPE exit code is recorded but does not discard physically valid
output, and a valid ACPYPE result is never overwritten.

Exit codes
----------
0  success
2  usage / environment error (bad input, missing converter/interpreter)
3  ACPYPE failed and left no usable Amber files
4  ACPYPE reported success but produced no usable files
5  charge validation failed: |sum(q) - net_charge| > 1e-4 e
6  GMX files are internally inconsistent (ITP/GRO atom counts differ)
7  ParmEd conversion failed or produced unusable files
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

CHARGE_TOL = 1e-4
EXIT_USAGE = 2
EXIT_ACPYPE_NO_AMBER = 3
EXIT_NO_OUTPUT = 4
EXIT_CHARGE = 5
EXIT_INCONSISTENT = 6
EXIT_CONVERTER = 7


# ---------------------------------------------------------------------------
# Argument inspection (we must pass every argument through to ACPYPE unchanged)
# ---------------------------------------------------------------------------
def _option(argv: Sequence[str], *names: str) -> Optional[str]:
    """Return the value of the first matching option (short or long form)."""
    for index, arg in enumerate(argv):
        if arg in names:
            if index + 1 < len(argv):
                return argv[index + 1]
            return None
        for name in names:
            if name.startswith("--") and arg.startswith(name + "="):
                return arg.split("=", 1)[1]
    return None


def _as_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Environment / subprocess helpers
# ---------------------------------------------------------------------------
def _amber_home() -> Optional[Path]:
    """Return the AmberTools prefix that actually contains antechamber."""
    candidates: List[Path] = []
    env_home = os.environ.get("AMBERHOME")
    if env_home:
        candidates.append(Path(env_home))
    candidates.append(Path(sys.executable).resolve().parents[1])
    candidates.append(Path("/opt/conda/envs/amber"))
    for candidate in candidates:
        if (candidate / "bin" / "antechamber").is_file():
            return candidate
    return candidates[0] if candidates else None


def _amber_env() -> dict:
    env = os.environ.copy()
    home = _amber_home()
    if home is not None:
        env["AMBERHOME"] = str(home)
        env["PATH"] = str(home / "bin") + os.pathsep + env.get("PATH", "")
        existing_lib = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = str(home / "lib") + (
            os.pathsep + existing_lib if existing_lib else ""
        )
    return env


def _run_acpype(argv: Sequence[str]) -> int:
    command = [sys.executable, "-m", "acpype.cli", *argv]
    print("[acpype-wrapper] running: " + " ".join(command), flush=True)
    try:
        return subprocess.run(command, env=_amber_env(), check=False).returncode
    except OSError as exc:
        print(f"[acpype-wrapper] ERROR: could not launch ACPYPE: {exc}", file=sys.stderr, flush=True)
        return EXIT_USAGE


def _converter_path() -> Optional[Path]:
    candidates = [
        Path("/usr/local/bin/convert_amber_to_gmx.py"),
        Path(__file__).resolve().with_name("convert_amber_to_gmx.py"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def _base_python() -> str:
    candidate = Path("/opt/conda/bin/python")
    return str(candidate) if candidate.is_file() else sys.executable


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------
def _itp_atoms(itp_path: Path) -> Tuple[int, float]:
    """Return (atom count, charge sum) from the ``[ atoms ]`` block."""
    count = 0
    total = 0.0
    in_atoms = False
    for raw in itp_path.read_text(encoding="utf-8", errors="replace").splitlines():
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
            int(parts[0])
            charge = float(parts[6])
        except ValueError:
            continue
        count += 1
        total += charge
    return count, total


def _gro_atoms(gro_path: Path) -> int:
    """Return the atom count from the second line of a GROMACS .gro file."""
    lines = gro_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        return 0
    try:
        return int(lines[1].strip())
    except ValueError:
        return 0


def _acpype_dir(root: str) -> Path:
    return Path(f"{root}.acpype")


def _gmx_state(root: str) -> Tuple[bool, int, int, float, str]:
    """Return (valid, itp_atoms, gro_atoms, charge_sum, detail)."""
    directory = _acpype_dir(root)
    itp = directory / f"{root}_GMX.itp"
    gro = directory / f"{root}_GMX.gro"
    if not itp.is_file() or not gro.is_file():
        return False, 0, 0, 0.0, "missing ITP/GRO"
    if itp.stat().st_size == 0 or gro.stat().st_size == 0:
        return False, 0, 0, 0.0, "empty ITP/GRO"
    n_itp, charge_sum = _itp_atoms(itp)
    n_gro = _gro_atoms(gro)
    if n_itp <= 0:
        return False, 0, n_gro, charge_sum, "no [ atoms ] entries in ITP"
    if n_gro <= 0:
        return False, n_itp, 0, charge_sum, "invalid GRO atom count"
    if n_itp != n_gro:
        return False, n_itp, n_gro, charge_sum, f"ITP/GRO atom mismatch ({n_itp} != {n_gro})"
    return True, n_itp, n_gro, charge_sum, f"atoms={n_itp} charge_sum={charge_sum:+.6f}"


def _amber_files(root: str) -> Optional[Tuple[Path, Path]]:
    """Return (prmtop, inpcrd) only when both look complete."""
    directory = _acpype_dir(root)
    prmtop = directory / f"{root}_AC.prmtop"
    inpcrd = directory / f"{root}_AC.inpcrd"
    if not prmtop.is_file() or not inpcrd.is_file():
        return None
    if prmtop.stat().st_size == 0 or inpcrd.stat().st_size == 0:
        return None
    text = prmtop.read_text(encoding="utf-8", errors="replace")
    required = ("%FLAG POINTERS", "%FLAG ATOM_NAME", "%FLAG CHARGE")
    if not all(flag in text for flag in required):
        return None
    if len(inpcrd.read_text(encoding="utf-8", errors="replace").splitlines()) < 2:
        return None
    return prmtop, inpcrd


def _run_converter(root: str, net_charge: Optional[int]) -> int:
    converter = _converter_path()
    if converter is None:
        print("[acpype-wrapper] ERROR: convert_amber_to_gmx.py not found", file=sys.stderr, flush=True)
        return EXIT_USAGE
    command = [_base_python(), str(converter), root, "--dir", str(_acpype_dir(root))]
    if net_charge is not None:
        command.extend(["--net-charge", str(int(net_charge))])
    print("[acpype-wrapper] running: " + " ".join(command), flush=True)
    try:
        return subprocess.run(command, check=False).returncode
    except OSError as exc:
        print(f"[acpype-wrapper] ERROR: could not launch converter: {exc}", file=sys.stderr, flush=True)
        return EXIT_USAGE


def _charge_target(charge_sum: float, net_charge: Optional[int]) -> int:
    # ACPYPE itself snaps to int(round(total)); mirror that when -n is absent.
    return int(net_charge) if net_charge is not None else int(round(charge_sum))


def _charge_ok(charge_sum: float, target: int) -> bool:
    return abs(charge_sum - float(target)) <= CHARGE_TOL


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv: Sequence[str]) -> int:
    resname = _option(argv, "-b", "--basename")
    input_path = _option(argv, "-i", "--input")
    net_charge = _as_int(_option(argv, "-n", "--net_charge"))
    outtop = (_option(argv, "-o", "--outtop") or "all").lower()

    if resname is None:
        if input_path is None:
            # -p/-x (amb2gmx) mode: output naming is not ours to infer.
            print("[acpype-wrapper] no -i/-b; passing through without fallback", flush=True)
            return _run_acpype(argv)
        # ACPYPE derives the root from the input stem when -b is absent.
        resname = Path(input_path).stem

    print(
        f"[acpype-wrapper] root={resname} input={input_path or '-'} "
        f"net_charge={net_charge if net_charge is not None else 'auto'} outtop={outtop}",
        flush=True,
    )

    if input_path is not None and not Path(input_path).is_file():
        print(f"[acpype-wrapper] ERROR: input file not found: {input_path}", file=sys.stderr, flush=True)
        return EXIT_USAGE

    acpype_rc = _run_acpype(argv)
    print(f"[acpype-wrapper] acpype_rc={acpype_rc}", flush=True)

    if outtop not in ("gmx", "all"):
        # No GROMACS pair is requested; propagate ACPYPE's result untouched.
        print(f"[acpype-wrapper] outtop={outtop}: no GROMACS fallback", flush=True)
        return acpype_rc

    amber = _amber_files(resname)
    gmx_valid, n_itp, n_gro, charge_sum, detail = _gmx_state(resname)
    print(
        f"[acpype-wrapper] state: amber={'yes' if amber else 'no'} "
        f"gmx={'valid' if gmx_valid else 'invalid'} ({detail})",
        flush=True,
    )

    target = _charge_target(charge_sum, net_charge)
    if gmx_valid and _charge_ok(charge_sum, target):
        print(
            f"[acpype-wrapper] accepted ACPYPE output unchanged; "
            f"charge {charge_sum:+.6f} e vs target {target:+d} e (rc={acpype_rc})",
            flush=True,
        )
        return 0
    if gmx_valid:
        print(
            f"[acpype-wrapper] ACPYPE output is not neutral: {charge_sum:+.6f} e "
            f"vs target {target:+d} e; re-balancing from Amber files",
            file=sys.stderr,
            flush=True,
        )

    if amber is None:
        if acpype_rc != 0:
            print(
                f"[acpype-wrapper] ERROR: ACPYPE failed (rc={acpype_rc}) and left no complete "
                "Amber files to convert",
                file=sys.stderr,
                flush=True,
            )
            return EXIT_ACPYPE_NO_AMBER
        if gmx_valid:
            target = _charge_target(charge_sum, net_charge)
            print(
                f"[acpype-wrapper] ERROR: charge validation failed: {charge_sum:+.6f} e "
                f"vs target {target:+d} e, and no Amber files exist to re-balance",
                file=sys.stderr,
                flush=True,
            )
            return EXIT_CHARGE
        print(
            "[acpype-wrapper] ERROR: ACPYPE reported success but produced neither a valid "
            "GMX pair nor complete Amber files",
            file=sys.stderr,
            flush=True,
        )
        return EXIT_NO_OUTPUT

    # Conversion path: ACPYPE failed after writing Amber files, or its GMX
    # output was unusable / non-neutral.
    converter_rc = _run_converter(resname, net_charge)
    if converter_rc != 0:
        print(
            f"[acpype-wrapper] ERROR: ParmEd conversion failed (rc={converter_rc})",
            file=sys.stderr,
            flush=True,
        )
        return EXIT_CONVERTER

    gmx_valid, n_itp, n_gro, charge_sum, detail = _gmx_state(resname)
    if not gmx_valid:
        print(
            f"[acpype-wrapper] ERROR: converted files are unusable ({detail})",
            file=sys.stderr,
            flush=True,
        )
        return EXIT_INCONSISTENT

    target = _charge_target(charge_sum, net_charge)
    if not _charge_ok(charge_sum, target):
        print(
            f"[acpype-wrapper] ERROR: charge validation failed after conversion: "
            f"{charge_sum:+.6f} e vs target {target:+d} e (tolerance {CHARGE_TOL:g} e)",
            file=sys.stderr,
            flush=True,
        )
        return EXIT_CHARGE

    print(
        f"[acpype-wrapper] SUCCESS: validated GMX output; atoms={n_itp} "
        f"charge={charge_sum:+.6f} e target={target:+d} e",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
