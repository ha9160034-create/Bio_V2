#!/usr/bin/env python3
"""Write or refresh ``md_manifest.json`` for a finished GROMACS MD run.

This script is intentionally standalone (standard library only) so it can be
shipped inside the exported reproducibility bundle and executed next to
``run_md.sh`` / ``run_md.bat``.  It lets the Streamlit app discover runs that
were executed outside the in-app orchestrator.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

PRODUCTION_FILES = ("prod_md.tpr", "prod_md.xtc")
REPLICATE_PATTERN = re.compile(r"^rep(\d+)_prod_md\.tpr$")


def count_replicates(run_dir: Path) -> int:
    """Count the primary plus every ``repN_prod_md.tpr`` production replica."""
    indices = {
        int(match.group(1))
        for path in run_dir.glob("rep*_prod_md.tpr")
        if (match := REPLICATE_PATTERN.match(path.name))
    }
    return 1 + len(indices)


def parse_mdp(path: Path) -> dict[str, str]:
    params: dict[str, str] = {}
    if not path.is_file():
        return params
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        params[key.strip()] = value.strip()
    return params


def detect_gromacs_version(run_dir: Path) -> str:
    for name in ("md.log", "prod_md.log", "mdrun_prod_md.log", "npt.log"):
        candidate = run_dir / name
        if not candidate.is_file():
            continue
        text = candidate.read_text(encoding="utf-8", errors="replace")
        match = re.search(r"GROMACS version:\s*(\S+)", text)
        if match:
            return match.group(1)
    return "unknown"


def production_finished(run_dir: Path) -> bool:
    """Require a normal mdrun termination, not merely the presence of files.

    A killed production job leaves ``prod_md.tpr`` and a partial ``prod_md.xtc``
    behind.  Marking that run as completed would let the UI load a truncated
    trajectory, so the finalizer also requires GROMACS' completion banner.
    """
    log_path = run_dir / "prod_md.log"
    if not log_path.is_file():
        return False
    text = log_path.read_text(encoding="utf-8", errors="replace")
    return "Finished mdrun" in text


def infer_ligand_resname(run_dir: Path) -> str | None:
    """Read the ligand residue name from an ACPYPE ``*_GMX.gro`` file.

    External runs have no initial manifest, so the saved-run loader would
    otherwise fall back to ``LIG`` and fail to highlight the real ligand
    (for example ``ERLOT`` for erlotinib).
    """
    for path in sorted(run_dir.glob("*_GMX.gro")):
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        if len(lines) < 3:
            continue
        resname = lines[2][5:10].strip()
        if resname:
            return resname
    return None


def finalize(run_dir: Path) -> dict:
    manifest_path = run_dir / "md_manifest.json"
    manifest: dict = {}
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            manifest = {}

    completed = all((run_dir / name).is_file() for name in PRODUCTION_FILES) and production_finished(run_dir)
    mdp = parse_mdp(run_dir / "md.mdp")
    try:
        nsteps = int(float(mdp.get("nsteps", "0")))
        dt = float(mdp.get("dt", "0.002"))
        duration_ns = nsteps * dt / 1000.0 if nsteps > 0 else None
    except ValueError:
        dt = 0.002
        duration_ns = None

    parameters = manifest.setdefault("parameters", {})
    if duration_ns is not None:
        parameters["production_duration_ns"] = duration_ns
    parameters["timestep_ps"] = dt
    parameters["n_replicates"] = count_replicates(run_dir)

    engine = manifest.setdefault("engine", {})
    engine.setdefault("name", "GROMACS")
    engine.setdefault("version", detect_gromacs_version(run_dir))

    topology = manifest.setdefault("topology", {})
    if not topology.get("ligand_resname"):
        inferred = infer_ligand_resname(run_dir)
        if inferred:
            topology["ligand_resname"] = inferred
    if not topology.get("ligand_itp"):
        for path in sorted(run_dir.glob("*_GMX.itp")):
            topology["ligand_itp"] = path.name
            break

    manifest.setdefault("run_id", run_dir.name)
    manifest.setdefault("protocol_reference", "wiki/gromacs-protein-ligand-protocol.md")
    manifest["status"] = "completed" if completed else "failed"
    if completed:
        manifest.pop("error", None)
    else:
        manifest.setdefault(
            "error",
            "Production did not finish: prod_md.tpr/prod_md.xtc are missing or "
            "prod_md.log has no 'Finished mdrun' banner.",
        )
    manifest["artifacts"] = sorted(path.name for path in run_dir.glob("*.xvg"))
    manifest["finalized_by"] = "finalize_md_run.py"

    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="MD run directory to finalize")
    args = parser.parse_args(argv)

    run_dir = Path(args.run_dir).resolve()
    if not run_dir.is_dir():
        print(f"Run directory not found: {run_dir}")
        return 2

    manifest = finalize(run_dir)
    print(f"Wrote {run_dir / 'md_manifest.json'} (status={manifest['status']})")
    return 0 if manifest["status"] == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
