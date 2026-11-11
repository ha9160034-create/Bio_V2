"""verify_stack.py — fail-closed verification of the BIOBIO runtime stack.

Designed from:
  * PLAN.md §1.1 / §1.2 / Phase 1 acceptance criteria
  * wiki/gromacs-protein-ligand-protocol.md (GROMACS + amber99sb-ildn defaults)
  * wiki/force-fields.md (GAFF2 / AMBER compatibility; fudgeLJ 0.5, fudgeQQ 0.8333)

Exit codes:
  0  every check passed
  1  at least one check failed (the build must be rejected)

Usage:
  python scripts/verify_stack.py [--expected /path/to/expected_versions.json]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

AMBER_PYTHON = Path("/opt/conda/envs/amber/bin/python")
GMX = shutil.which("gmx") or "/usr/local/bin/gmx"


def run(cmd: list[str], timeout: int = 120) -> tuple[int, str]:
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        return proc.returncode, (proc.stdout or "") + (proc.stderr or "")
    except Exception as exc:  # pragma: no cover - environment failure
        return -1, f"{type(exc).__name__}: {exc}"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_defaults(forcefield_itp: Path) -> tuple[int, float, float] | None:
    """Return (comb_rule, fudgeLJ, fudgeQQ) from a GROMACS forcefield.itp."""
    lines = forcefield_itp.read_text(errors="replace").splitlines()
    in_defaults = False
    for line in lines:
        stripped = line.strip()
        if stripped.lower().startswith("[ defaults ]"):
            in_defaults = True
            continue
        if in_defaults:
            if stripped.startswith(";") or not stripped:
                continue
            if stripped.startswith("["):
                return None
            parts = stripped.split()
            if len(parts) < 5:
                return None
            return int(parts[1]), float(parts[3]), float(parts[4])
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--expected",
        default=str(Path(__file__).resolve().parent.parent / "expected_versions.json"),
    )
    args = parser.parse_args()

    expected = json.loads(Path(args.expected).read_text(encoding="utf-8"))
    versions = expected["versions"]

    checks: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str) -> None:
        checks.append((name, bool(ok), detail))

    # --- Python ---
    py = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    check("python", py == versions["python"], f"{py} (expected {versions['python']})")

    # --- RDKit ---
    try:
        import rdkit

        check("rdkit", rdkit.__version__ == versions["rdkit"], f"{rdkit.__version__}")
    except Exception as exc:
        check("rdkit", False, f"import failed: {exc}")

    # --- OpenBabel (wheel) ---
    try:
        from openbabel import openbabel as ob

        ob_version = ob.OBReleaseVersion()
        check("openbabel(py)", ob_version == versions["openbabel"], ob_version)
    except Exception as exc:
        check("openbabel(py)", False, f"import failed: {exc}")

    obabel = shutil.which("obabel")
    if obabel:
        rc, out = run([obabel, "-V"])
        m = re.search(r"Open Babel\s+([0-9.]+)", out)
        check(
            "openbabel(cli)",
            rc == 0 and bool(m) and m.group(1) == versions["openbabel"],
            (m.group(1) if m else out.strip()[:80]),
        )
    else:
        check("openbabel(cli)", False, "obabel not found in PATH")

    # --- ParmEd / NumPy / SciPy / Streamlit ---
    for module, key in (
        ("parmed", "parmed"),
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("streamlit", "streamlit"),
    ):
        try:
            mod = __import__(module)
            got = getattr(mod, "__version__", "unknown")
            check(module, got == versions[key], got)
        except Exception as exc:
            check(module, False, f"import failed: {exc}")

    # --- GROMACS ---
    rc, out = run([GMX, "--version"])
    gmx_version = None
    simd = None
    for line in out.splitlines():
        if "GROMACS version:" in line:
            gmx_version = line.split(":", 1)[1].strip().split("-")[0]
        if "SIMD instructions:" in line:
            simd = line.split(":", 1)[1].strip()
    check(
        "gromacs",
        rc == 0 and gmx_version == versions["gromacs"],
        f"{gmx_version} (expected {versions['gromacs']})",
    )
    check(
        "gromacs.simd",
        simd in expected["gromacs_simd_allowed"],
        f"{simd} (allowed {expected['gromacs_simd_allowed']})",
    )

    # --- Vina: must be the pinned binary, and must not be shadowed by apt ---
    vina = Path(expected["vina"]["path"])
    if vina.is_file():
        digest = sha256(vina)
        rc, out = run([str(vina), "--version"])
        m = re.search(r"AutoDock Vina v([0-9.]+)", out)
        got_version = f"AutoDock Vina v{m.group(1)}" if m else out.strip()[:80]
        check(
            "vina.version",
            rc == 0 and got_version == expected["vina"]["version"],
            f"{got_version} (expected {expected['vina']['version']})",
        )
        check(
            "vina.sha256",
            digest == expected["vina"]["sha256"],
            f"{digest[:16]}… (expected {expected['vina']['sha256'][:16]}…)",
        )
    else:
        check("vina", False, f"missing {vina}")

    for forbidden in expected["forbidden_paths"]:
        check(
            f"no duplicate {forbidden}",
            not Path(forbidden).exists(),
            "absent" if not Path(forbidden).exists() else "PRESENT — duplicate Vina",
        )
    check(
        "vina not on PATH",
        shutil.which("vina") is None,
        f"{shutil.which('vina')}",
    )

    # --- AmberTools / ACPYPE (amber env) ---
    antechamber = Path("/opt/conda/envs/amber/bin/antechamber")
    if antechamber.is_file():
        rc, out = run([str(antechamber), "-h"])
        check("ambertools.antechamber", "Usage: antechamber" in out, "responds to -h")
    else:
        check("ambertools.antechamber", False, "missing binary")

    if AMBER_PYTHON.is_file():
        rc, out = run(
            [str(AMBER_PYTHON), "-c", "import acpype; print(acpype.__version__)"]
        )
        got = out.strip().splitlines()[-1] if out.strip() else ""
        check("acpype", rc == 0 and got == versions["acpype"], got)
    else:
        check("acpype", False, f"missing {AMBER_PYTHON}")

    # --- Force-field defaults (justifies removing [ defaults ] from ligand ITP) ---
    ff_candidates = [
        Path("/opt/conda/share/gromacs/top") / f"{expected['forcefield']['name']}.ff" / "forcefield.itp",
    ]
    gmx_data = os.environ.get("GMXDATA")
    if gmx_data:
        ff_candidates.insert(0, Path(gmx_data) / "top" / f"{expected['forcefield']['name']}.ff" / "forcefield.itp")
    ff_file = next((p for p in ff_candidates if p.is_file()), None)
    if ff_file:
        defaults = parse_defaults(ff_file)
        ff = expected["forcefield"]
        ok = (
            defaults is not None
            and defaults[0] == ff["comb_rule"]
            and abs(defaults[1] - ff["fudgeLJ"]) < 1e-9
            and abs(defaults[2] - ff["fudgeQQ"]) < 1e-9
        )
        check(
            "amber99sb-ildn defaults",
            ok,
            f"{defaults} (expected ({ff['comb_rule']}, {ff['fudgeLJ']}, {ff['fudgeQQ']}))",
        )
    else:
        check("amber99sb-ildn defaults", False, "forcefield.itp not found")

    # --- Required application files (guards the rm -rf regression, issue #10) ---
    for rel in expected["required_files"]:
        path = Path(rel)
        check(f"file {rel}", path.is_file(), "present" if path.is_file() else "MISSING")

    # --- Report ---
    width = max(len(name) for name, _, _ in checks) + 2
    failures = 0
    print("BIOBIO stack verification")
    print("=" * (width + 60))
    for name, ok, detail in checks:
        status = "PASS" if ok else "FAIL"
        if not ok:
            failures += 1
        print(f"{status:4}  {name:<{width}}  {detail}")
    print("=" * (width + 60))
    print(f"{len(checks) - failures}/{len(checks)} checks passed")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
