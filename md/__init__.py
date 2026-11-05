"""
Molecular Dynamics (MD) Simulation & Analysis Package
=====================================================

A modular biophysical platform for all-atom molecular dynamics simulation of
protein-ligand complexes derived from molecular docking runs.

Scientific Foundations & Wiki Traceability:
-------------------------------------------
This package strictly adheres to the protocols and biophysical guardrails defined in:
- `wiki/gromacs-protein-ligand-protocol.md`: GROMACS Protein-Ligand Simulation Protocol (Justin Lemkul / LiveCoMS).
- `wiki/open-force-field-and-ligand-parameterization.md`: Small-Molecule Force Fields (OpenFF Sage, CGenFF, ACPYPE/GAFF).
- `wiki/simulation-and-benchmark-best-practices.md`: LiveCoMS Foundations & Benchmark Best Practices.
- `wiki/mm-pbsa-and-mm-gbsa.md`: Endpoint Binding Free Energy Rescoring & 1-Trajectory Approach.
- `wiki/force-fields.md`: Force Field Architectures & Nonbonded Formulations.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MD_RUNS_DIR = PROJECT_ROOT / "md_runs"

__version__ = "1.0.0"
__all__ = [
    "PROJECT_ROOT",
    "MD_RUNS_DIR",
]

