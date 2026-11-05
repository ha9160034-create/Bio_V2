"""
Molecular Dynamics Simulation Pipeline & MDP Configuration Layer
================================================================

Scientific Traceability & Wiki Documentation:
---------------------------------------------
- `wiki/gromacs-protein-ligand-protocol.md`: Lemkul's Step-by-Step Simulation Protocol
  (Box definition, solvation, neutralization, EM, NVT, NPT, Production MD).
- `wiki/simulation-and-benchmark-best-practices.md`: LiveCoMS Foundations in Molecular Simulations
  (Integrator limits, thermostat audit, barostat ensemble rules, replicate sampling).
- `wiki/open-force-field-and-ligand-parameterization.md`: Force field compatibility & cutoff rules.
- `wiki/mm-pbsa-and-mm-gbsa.md`: Production trajectory stability prerequisites.

Biophysical Guardrails (per AGENTS.md):
---------------------------------------
1. Unit Cell: Default to rhombic dodecahedral geometry (`-bt dodecahedron -d 1.0`) to save ~23% solvent volume.
2. Restraints & Thermostats: Restrain ligand heavy atoms during solvent equilibration, and merge
   protein and ligand into a single thermostat group (`tc-grps = Protein_Ligand Water_and_ions`).
3. Ensemble Validity: Use Bussi velocity-rescaling (`v-rescale`) or Nosé–Hoover. NEVER use the
   Berendsen thermostat or barostat in production MD simulations.
4. Trajectory Processing: Always apply periodic boundary condition (PBC) re-wrapping
   (`trjconv -pbc mol -ur compact` and `-fit rot+trans`) before quantitative analysis.

No-fabrication policy
---------------------
Earlier versions silently synthesised a "trajectory" when GROMACS was absent.
That behaviour has been removed. If GROMACS cannot run, the pipeline raises
:class:`SimulationError` and writes a failed manifest. Only files produced by
GROMACS are reported as results.
"""

from __future__ import annotations

import os
import re
import math
import json
import uuid
import shutil
import hashlib
import datetime
import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any, Callable

from md.engine import (
    GromacsBackend,
    GromacsExecutionError,
    detect_gromacs,
    generate_standalone_scripts,
    run_gmx_command,
    solvent_coordinate_file,
)
from md.indexing import parse_gro_atoms, write_index_from_gro
from md.topology import TopologyResult
from md.analysis import parse_xvg_file, summarize_timeseries

logger = logging.getLogger("md.simulation")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MD_RUNS_DIR = PROJECT_ROOT / "md_runs"


class SimulationError(RuntimeError):
    """Raised when a real GROMACS simulation cannot be completed."""


@dataclass
class SimulationConfig:
    run_id: str = ""
    prefix: str = "md"
    temperature_k: float = 300.0
    pressure_bar: float = 1.0
    timestep_ps: float = 0.002       # 2.0 fs
    em_max_steps: int = 50000
    em_tol: float = 1000.0           # Fmax < 1000 kJ/(mol*nm)
    nve_duration_ps: float = 10.0    # short NVE drift check before NVT
    nvt_duration_ps: float = 100.0   # 100 ps
    npt_duration_ps: float = 100.0   # 100 ps
    prod_duration_ns: float = 1.0    # 1.0 ns default
    box_type: str = "dodecahedron"   # Rhombic Dodecahedron (saves ~23% water volume)
    box_distance_nm: float = 1.0     # 1.0 nm minimum distance to wall
    water_model: str = "tip3p"
    salt_concentration_molar: float = 0.15  # 0.15 M physiological NaCl
    force_field: str = "amber99sb-ildn"
    ligand_restraint_fc: float = 1000.0  # kJ/(mol*nm^2)
    random_seed: int = 42
    grompp_maxwarn: int = 2
    # Number of independent production trajectories.  With >1, extra
    # trajectories start from the same equilibrated NPT checkpoint but with
    # independent Maxwell-Boltzmann velocities (see manifest replication_note).
    n_replicates: int = 1
    # Fraction of the production trajectory discarded as equilibration/burn-in
    # when reporting mean RMSD/H-bond metrics.
    analysis_burn_in_fraction: float = 0.2


def generate_mdp_ions() -> str:
    """Generate MDP parameters for ion addition and net charge neutralization."""
    return """; MDP for ion addition (genion)
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
cutoff-scheme = Verlet
nstlist     = 10
rlist       = 1.0
coulombtype = cutoff
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
"""


def _nonbonded_section(force_field: str) -> str:
    if "charmm" in force_field.lower():
        return """vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw            = 1.2
coulombtype     = PME
rcoulomb        = 1.2
pme_order       = 4
fourierspacing  = 0.16
DispCorr        = no"""
    return """vdwtype         = cutoff
vdw-modifier    = none
rvdw            = 1.0
coulombtype     = PME
rcoulomb        = 1.0
pme_order       = 4
fourierspacing  = 0.16
DispCorr        = EnerPres"""


def generate_mdp_minimization(force_field: str = "amber99sb-ildn", em_tol: float = 1000.0, nsteps: int = 50000) -> str:
    """
    Generate MDP parameters for Energy Minimization using Steepest Descents.
    Implements force-switch nonbonded settings for CHARMM36 and AMBER conventions.
    """
    return f"""; ====================================================================
; Energy Minimization (Steepest Descents)
; Protocol: Lemkul & LiveCoMS Foundations (wiki/gromacs-protein-ligand-protocol.md)
; ====================================================================
integrator      = steep
emtol           = {em_tol:.1f}
emstep          = 0.01
nsteps          = {nsteps:d}
nstlist         = 20
cutoff-scheme   = Verlet
rlist           = 1.0
{_nonbonded_section(force_field)}
pbc             = xyz
"""


def generate_mdp_nvt(
    temperature_k: float = 300.0,
    duration_ps: float = 100.0,
    timestep_ps: float = 0.002,
    force_field: str = "amber99sb-ildn",
    random_seed: int = 42,
) -> str:
    """
    Generate MDP parameters for NVT thermalization (canonical ensemble).
    Guardrails:
    - Bussi velocity-rescaling thermostat (tcoupl = v-rescale)
    - Merged thermostat group: tc-grps = Protein_Ligand Water_and_ions
    - Ligand heavy atoms & protein position restrained (define = -DPOSRES -DPOSRES_LIG)
    """
    nsteps = int(round(duration_ps / timestep_ps))
    return f"""; ====================================================================
; NVT Thermalization ({duration_ps:.0f} ps, Canonical Ensemble)
; Protocol: Lemkul & LiveCoMS Foundations (wiki/gromacs-protein-ligand-protocol.md)
; ====================================================================
title                   = NVT Equilibration with Ligand Heavy-Atom Restraints
define                  = -DPOSRES -DPOSRES_LIG

; Run controls
integrator              = md
dt                      = {timestep_ps:.3f}
nsteps                  = {nsteps:d}

; Output control
nstxout-compressed      = 1000
nstenergy               = 1000
nstlog                  = 1000

; Bond parameters (LINCS for a 2.0 fs timestep per LiveCoMS)
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

; Nonbonded cutoffs
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
{_nonbonded_section(force_field)}

; Temperature coupling (Bussi stochastic velocity-rescaling)
; CRITICAL GUARDRAIL: merged Protein_Ligand and Water_and_ions groups
tcoupl                  = v-rescale
tc-grps                 = Protein_Ligand Water_and_ions
tau_t                   = 0.1             0.1
ref_t                   = {temperature_k:.1f}           {temperature_k:.1f}

; Pressure coupling is OFF in NVT
pcoupl                  = no

; Periodic boundary conditions
pbc                     = xyz

; Velocity generation
gen_vel                 = yes
gen_temp                = {temperature_k:.1f}
gen_seed                = {random_seed:d}
"""


def generate_mdp_nve(
    duration_ps: float = 10.0,
    temperature_k: float = 300.0,
    timestep_ps: float = 0.002,
    force_field: str = "amber99sb-ildn",
    random_seed: int = 42,
) -> str:
    """Short NVE run to check integrator/cutoff energy drift (LiveCoMS).

    No thermostat and no barostat are active.  A short NVE segment is the
    standard way to verify that the chosen timestep and nonbonded settings do
    not introduce systematic total-energy drift before production.
    """
    nsteps = max(10, int(round(duration_ps / timestep_ps)))
    return f"""; ====================================================================
; NVE Drift Check ({duration_ps:.1f} ps, microcanonical)
; Protocol: wiki/simulation-and-benchmark-best-practices.md
; ====================================================================
title                   = NVE drift check (no thermostat, no barostat)
integrator              = md
dt                      = {timestep_ps:.3f}
nsteps                  = {nsteps:d}
continuation            = yes

; No temperature or pressure coupling in NVE
tcoupl                  = no
pcoupl                  = no

; Bond parameters
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

; Nonbonded cutoffs (must match production)
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
{_nonbonded_section(force_field)}

; Periodic boundary conditions
pbc                     = xyz

; Output: energy every 10 steps, no compressed trajectory needed
nstenergy               = 10
nstlog                  = 100
nstxout-compressed      = 0

; Initial velocities: continue from the equilibrated NPT checkpoint
gen_vel                 = no
"""


def generate_mdp_npt(
    temperature_k: float = 300.0,
    pressure_bar: float = 1.0,
    duration_ps: float = 100.0,
    timestep_ps: float = 0.002,
    force_field: str = "amber99sb-ildn",
) -> str:
    """
    Generate MDP parameters for NPT isobaric-isothermal equilibration.
    Berendsen is permitted strictly during equilibration; production uses
    Parrinello-Rahman (see generate_mdp_production).
    """
    nsteps = int(round(duration_ps / timestep_ps))
    return f"""; ====================================================================
; NPT Equilibration ({duration_ps:.0f} ps, Isobaric-Isothermal Ensemble)
; Protocol: Lemkul & LiveCoMS Foundations (wiki/gromacs-protein-ligand-protocol.md)
; ====================================================================
title                   = NPT Equilibration with Ligand Heavy-Atom Restraints
define                  = -DPOSRES -DPOSRES_LIG

; Run controls
integrator              = md
dt                      = {timestep_ps:.3f}
nsteps                  = {nsteps:d}

; Output control
nstxout-compressed      = 1000
nstenergy               = 1000
nstlog                  = 1000

; Bond parameters
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

; Nonbonded cutoffs
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
{_nonbonded_section(force_field)}

; Temperature coupling (Bussi stochastic velocity-rescaling)
tcoupl                  = v-rescale
tc-grps                 = Protein_Ligand Water_and_ions
tau_t                   = 0.1             0.1
ref_t                   = {temperature_k:.1f}           {temperature_k:.1f}

; Pressure coupling (Berendsen allowed strictly during equilibration)
pcoupl                  = berendsen
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = {pressure_bar:.1f}
compressibility         = 4.5e-5
refcoord_scaling        = com

; Periodic boundary conditions
pbc                     = xyz

; Velocity generation
gen_vel                 = no
"""


def generate_mdp_production(
    duration_ns: float = 1.0,
    temperature_k: float = 300.0,
    pressure_bar: float = 1.0,
    timestep_ps: float = 0.002,
    force_field: str = "amber99sb-ildn",
    barostat: str = "Parrinello-Rahman",
    gen_vel: bool = False,
    random_seed: int = 42,
) -> str:
    """
    Generate MDP parameters for unrestrained Production Molecular Dynamics.
    CRITICAL BIOPHYSICAL GUARDRAIL (AGENTS.md & LiveCoMS):
    - Rejects Berendsen barostat/thermostat in production.
    - Uses Parrinello-Rahman barostat for a genuine isobaric ensemble.
    - All position restraints released.
    """
    if "berendsen" in barostat.lower():
        raise ValueError(
            "CRITICAL PROTOCOL VIOLATION: Berendsen barostat is strictly prohibited in production MD "
            "(see wiki/simulation-and-benchmark-best-practices.md & AGENTS.md). "
            "Use Parrinello-Rahman barostat."
        )

    duration_ps = duration_ns * 1000.0
    nsteps = int(round(duration_ps / timestep_ps))
    # Keep roughly 100 frames for a 1 ns run (10 ps spacing), never 0.
    nstxout = max(1, int(round(nsteps / 100.0)))
    continuation = "no" if gen_vel else "yes"
    velocity_block = (
        f"gen_vel                 = yes\ngen_temp                = {temperature_k:.1f}\n"
        f"gen_seed                = {random_seed:d}"
        if gen_vel
        else "gen_vel                 = no"
    )
    return f"""; ====================================================================
; Production Molecular Dynamics (Unrestrained, Parrinello-Rahman Barostat)
; Protocol: Lemkul & LiveCoMS Foundations (wiki/gromacs-protein-ligand-protocol.md)
; ====================================================================
title                   = Unrestrained Production MD ({duration_ns:.2f} ns)
; Restraints released

; Run controls
integrator              = md
dt                      = {timestep_ps:.3f}
nsteps                  = {nsteps:d}

; Output control
nstxout-compressed      = {nstxout}
compressed-x-grps       = System
nstenergy               = {nstxout}
nstlog                  = {nstxout}

; Bond parameters
continuation            = {continuation}
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4

; Nonbonded cutoffs
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.0
{_nonbonded_section(force_field)}

; Temperature coupling (Bussi stochastic velocity-rescaling)
tcoupl                  = v-rescale
tc-grps                 = Protein_Ligand Water_and_ions
tau_t                   = 0.1             0.1
ref_t                   = {temperature_k:.1f}           {temperature_k:.1f}

; Pressure coupling: PARRINELLO-RAHMAN IS MANDATORY FOR PRODUCTION
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = {pressure_bar:.1f}
compressibility         = 4.5e-5

; Periodic boundary conditions
pbc                     = xyz

; Velocity generation
{velocity_block}
"""


def compute_rhombic_dodecahedron_volume(
    box_distance_nm: float, solute_bounding_diameter_nm: float
) -> Tuple[float, float, float]:
    """
    Calculate unit cell volume and savings for a rhombic dodecahedron vs a cubic box.
    Returns (dodecahedron_vol_nm3, cubic_vol_nm3, savings_percentage).
    """
    d = solute_bounding_diameter_nm + (2.0 * box_distance_nm)
    cube_vol = d ** 3
    dodec_vol = 0.5 * math.sqrt(2.0) * (d ** 3)
    savings = ((cube_vol - dodec_vol) / cube_vol) * 100.0
    return dodec_vol, cube_vol, savings


def solute_bounding_diameter_nm(complex_gro: Path) -> Optional[float]:
    """Largest coordinate extent of the prepared complex, in nm.

    The dodecahedron volume depends on the actual solute diameter; using a
    hard-coded 5 nm value makes the reported box volume meaningless for small
    or large complexes.
    """
    try:
        atoms, _ = parse_gro_atoms(Path(complex_gro))
        if not atoms:
            return None
        xs = [atom[3] for atom in atoms]
        ys = [atom[4] for atom in atoms]
        zs = [atom[5] for atom in atoms]
        return max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
    except Exception:
        return None


def create_md_run_directory(prefix: str = "md") -> Tuple[str, Path]:
    """Create a timestamped MD run directory inside md_runs/."""
    MD_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S")
    uid = uuid.uuid4().hex[:6]
    run_id = f"{prefix}_{ts}_{uid}"
    run_dir = MD_RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_id, run_dir


def _software_versions(backend: Optional[GromacsBackend]) -> Dict[str, str]:
    """Record the software stack that produced the run."""
    import platform

    versions: Dict[str, str] = {
        "python": platform.python_version(),
        "gromacs": backend.version if backend else "not_run",
        "gromacs_backend": backend.kind if backend else "not_run",
    }
    try:
        import rdkit

        versions["rdkit"] = getattr(rdkit, "__version__", "unknown")
    except Exception:
        versions["rdkit"] = "unavailable"
    try:
        from openbabel import openbabel as ob

        versions["openbabel"] = ob.OBReleaseVersion()
    except Exception:
        versions["openbabel"] = "unavailable"
    return versions


def _sha256_file(path: Path) -> str:
    """SHA-256 of an MD artifact, for manifest provenance."""
    path = Path(path)
    if not path.is_file():
        return ""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _mean_std(values: List[float]) -> Tuple[Optional[float], Optional[float]]:
    if not values:
        return None, None
    n = len(values)
    mean = sum(values) / n
    if n < 2:
        return mean, 0.0
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    return mean, math.sqrt(var)


def parse_gmx_energy_xvg(path: Path) -> Dict[str, List[float]]:
    """Parse a multi-column ``gmx energy`` XVG file into {legend: values}."""
    path = Path(path)
    if not path.is_file():
        return {}
    legends: List[str] = []
    columns: Dict[str, List[float]] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("@ s") and "legend" in line:
                match = re.search(r'legend\s+"([^"]+)"', line)
                if match:
                    legends.append(match.group(1))
            elif not line or line.startswith(("#", "@")):
                continue
            else:
                parts = line.split()
                if len(parts) < 2:
                    continue
                values = []
                for token in parts[1:]:
                    try:
                        values.append(float(token))
                    except ValueError:
                        values.append(float("nan"))
                for i, val in enumerate(values):
                    if i < len(legends):
                        columns.setdefault(legends[i], []).append(val)
    return columns


def parse_hbond_num_xvg(path: Path) -> List[Tuple[float, float]]:
    """Parse ``gmx hbond -num`` into [(time_ps, hbond_count)]."""
    path = Path(path)
    if not path.is_file():
        return []
    points: List[Tuple[float, float]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    points.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    continue
    return points


class BiophysicalMDSimulator:
    """
    Real GROMACS MD executor for docked protein–ligand complexes.

    The class never fabricates coordinates, trajectories or energies. Every
    reported artifact must exist on disk and have been produced by GROMACS.
    """

    def __init__(
        self,
        config: SimulationConfig,
        output_dir: Path,
        backend: Optional[GromacsBackend] = None,
    ):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.backend = backend if backend is not None else detect_gromacs()

    # -- helpers ---------------------------------------------------------

    def _gmx(self, args: List[str], log_name: str, input_text: Optional[str] = None,
             timeout: Optional[int] = 1800) -> Tuple[int, str, str]:
        if self.backend is None:
            raise SimulationError(
                "GROMACS is not available (native or WSL). Simulation cannot start."
            )
        return run_gmx_command(
            args,
            cwd=self.output_dir,
            input_text=input_text,
            log_name=log_name,
            timeout=timeout,
            backend=self.backend,
        )

    def _mdrun(self, deffnm: str, stage_label: str) -> None:
        # Production MD may legitimately take hours: no timeout by default.
        self._gmx(["mdrun", "-v", "-deffnm", deffnm], log_name=f"mdrun_{deffnm}", timeout=None)

    def _stage_topology(self, topology: TopologyResult) -> TopologyResult:
        """Copy topology inputs into the run directory (self-contained run)."""
        out = self.output_dir

        def _copy(src: Path) -> Path:
            src = Path(src)
            dst = out / src.name
            if src.resolve() != dst.resolve():
                shutil.copy2(src, dst)
            return dst

        top = _copy(topology.topol_top)
        itp = _copy(topology.ligand_itp)
        posre = _copy(topology.ligand_posre)
        complex_gro = _copy(topology.complex_gro)
        complex_pdb = _copy(topology.complex_pdb)

        # pdb2gmx topol references the ligand ITP by basename; ensure it matches.
        text = top.read_text(encoding="utf-8", errors="replace")

        def _repl(match: "re.Match[str]") -> str:
            target = match.group(1)
            if ".ff/" in target or "forcefield.itp" in target or Path(target).name.startswith("posre"):
                return match.group(0)
            return f'#include "{itp.name}"'

        text = re.sub(r'#include\s+"([^"]+\.itp)"', _repl, text)
        text = text.replace("\r\n", "\n")
        top.write_text(text, encoding="utf-8", newline="\n")

        return TopologyResult(
            **{
                **topology.__dict__,
                "topol_top": top,
                "ligand_itp": itp,
                "ligand_posre": posre,
                "complex_gro": complex_gro,
                "complex_pdb": complex_pdb,
            }
        )

    def _write_mdp_files(self) -> Dict[str, Path]:
        cfg = self.config
        files = {
            "ions.mdp": generate_mdp_ions(),
            "em.mdp": generate_mdp_minimization(cfg.force_field, cfg.em_tol, cfg.em_max_steps),
            "nve.mdp": generate_mdp_nve(
                cfg.nve_duration_ps, cfg.temperature_k, cfg.timestep_ps,
                cfg.force_field, cfg.random_seed,
            ),
            "nvt.mdp": generate_mdp_nvt(
                cfg.temperature_k, cfg.nvt_duration_ps, cfg.timestep_ps,
                cfg.force_field, cfg.random_seed,
            ),
            "npt.mdp": generate_mdp_npt(
                cfg.temperature_k, cfg.pressure_bar, cfg.npt_duration_ps,
                cfg.timestep_ps, cfg.force_field,
            ),
            "md.mdp": generate_mdp_production(
                cfg.prod_duration_ns, cfg.temperature_k, cfg.pressure_bar,
                cfg.timestep_ps, cfg.force_field,
            ),
        }
        paths: Dict[str, Path] = {}
        for name, content in files.items():
            p = self.output_dir / name
            p.write_text(content, encoding="utf-8", newline="\n")
            paths[name] = p
        return paths

    # -- pipeline --------------------------------------------------------

    def run_pipeline(
        self,
        topology: TopologyResult,
        progress_callback: Optional[Callable[[str, float], None]] = None,
    ) -> Dict[str, Any]:
        """
        Execute the full Lemkul/LiveCoMS staged MD pipeline on a real topology:
        box setup, solvation, ions, index generation, EM, NVT, NPT, production,
        PBC post-processing, and quantitative trajectory analysis.
        """
        cfg = self.config
        out = self.output_dir
        n_replicates = max(1, min(int(cfg.n_replicates), 5))

        def update(stage: str, frac: float) -> None:
            if progress_callback:
                progress_callback(stage, frac)
            logger.info("MD Pipeline [%.1f%%]: %s", frac * 100.0, stage)

        if self.backend is None:
            manifest = self._write_failed_manifest(
                "GROMACS was not found natively or through WSL. No simulation was run.",
                stage="engine_detection",
            )
            raise SimulationError(manifest["error"])

        update("التحقق من ملفات الطوبولوجيا وإعداد دليل التشغيل...", 0.02)
        try:
            staged = self._stage_topology(topology)
            mdp = self._write_mdp_files()

            update("توليد سكربتات إعادة الإنتاج (Bash/Batch)...", 0.05)
            script_paths = generate_standalone_scripts(
                output_dir=out,
                complex_gro=staged.complex_gro.name,
                topol_top=staged.topol_top.name,
                index_ndx="index.ndx",
                box_type=cfg.box_type,
                box_distance=cfg.box_distance_nm,
                production_prefix="prod_md",
                ligand_resname=staged.ligand_resname,
                water_model=cfg.water_model,
                n_replicates=n_replicates,
                salt_concentration_molar=cfg.salt_concentration_molar,
            )
            # Ship the exact index generator used by the pipeline.
            shutil.copy2(Path(__file__).resolve().parent / "indexing.py", out / "make_index.py")
        except Exception as exc:
            self._write_failed_manifest(str(exc), stage="topology_staging")
            raise SimulationError(str(exc)) from exc

        warnings: List[str] = list(staged.warnings)
        nve_summary: Dict[str, Any] = {}
        replicate_stats: List[Dict[str, Any]] = []
        solute_diameter_nm = solute_bounding_diameter_nm(staged.complex_gro)
        if solute_diameter_nm is None or solute_diameter_nm <= 0:
            solute_diameter_nm = 5.0
            warnings.append(
                "Could not derive the solute diameter from complex.gro; "
                "box-volume savings are reported with a 5 nm fallback."
            )
        dodec_vol, cube_vol, savings_pct = compute_rhombic_dodecahedron_volume(
            cfg.box_distance_nm, solute_diameter_nm
        )
        failed_stage = "gromacs_execution"

        try:
            update("إنشاء صندوق اثنا عشري معيني (editconf)...", 0.10)
            self._gmx(
                ["editconf", "-f", staged.complex_gro.name, "-o", "newbox.gro",
                 "-bt", cfg.box_type, "-d", str(cfg.box_distance_nm)],
                log_name="editconf",
            )

            update("إذابة النظام في الماء الصريح (solvate)...", 0.18)
            self._gmx(
                ["solvate", "-cp", "newbox.gro", "-cs", solvent_coordinate_file(cfg.water_model),
                 "-p", staged.topol_top.name, "-o", "solv.gro"],
                log_name="solvate",
            )

            update("إضافة الأيونات الفسيولوجية (genion)...", 0.26)
            self._gmx(
                ["grompp", "-f", mdp["ions.mdp"].name, "-c", "solv.gro",
                 "-p", staged.topol_top.name, "-o", "ions.tpr",
                 "-maxwarn", str(cfg.grompp_maxwarn)],
                log_name="grompp_ions",
            )
            self._gmx(
                ["genion", "-s", "ions.tpr", "-o", "solv_ions.gro",
                 "-p", staged.topol_top.name, "-pname", "NA", "-nname", "CL",
                 "-neutral", "-conc", str(cfg.salt_concentration_molar)],
                log_name="genion",
                input_text="SOL\n",
            )

            update("بناء مجموعات الفهرس من النظام الفعلي (index.ndx)...", 0.32)
            failed_stage = "index_generation"
            group_sizes = write_index_from_gro(
                out / "solv_ions.gro", out / "index.ndx", ligand_resname=staged.ligand_resname
            )
            failed_stage = "gromacs_execution"
            logger.info("Index groups: %s", group_sizes)

            update("تقليل الطاقة (Energy Minimization)...", 0.40)
            self._gmx(
                ["grompp", "-f", mdp["em.mdp"].name, "-c", "solv_ions.gro",
                 "-p", staged.topol_top.name, "-o", "em.tpr",
                 "-maxwarn", str(cfg.grompp_maxwarn)],
                log_name="grompp_em",
            )
            self._mdrun("em", "EM")

            update("التوازن الحراري NVT (v-rescale)...", 0.55)
            self._gmx(
                ["grompp", "-f", mdp["nvt.mdp"].name, "-c", "em.gro", "-r", "em.gro",
                 "-p", staged.topol_top.name, "-n", "index.ndx", "-o", "nvt.tpr",
                 "-maxwarn", str(cfg.grompp_maxwarn)],
                log_name="grompp_nvt",
            )
            self._mdrun("nvt", "NVT")

            update("التوازن الضغطي NPT (Berendsen أثناء التوازن فقط)...", 0.70)
            self._gmx(
                ["grompp", "-f", mdp["npt.mdp"].name, "-c", "nvt.gro", "-t", "nvt.cpt",
                 "-r", "nvt.gro", "-p", staged.topol_top.name, "-n", "index.ndx",
                 "-o", "npt.tpr", "-maxwarn", str(cfg.grompp_maxwarn)],
                log_name="grompp_npt",
            )
            self._mdrun("npt", "NPT")

            if cfg.nve_duration_ps > 0:
                update("فحص انحراف الطاقة NVE بعد التوازن (بدون ترموستات/باروستات)...", 0.78)
                self._gmx(
                    ["grompp", "-f", mdp["nve.mdp"].name, "-c", "npt.gro", "-t", "npt.cpt",
                     "-p", staged.topol_top.name, "-n", "index.ndx", "-o", "nve.tpr",
                     "-maxwarn", str(cfg.grompp_maxwarn)],
                    log_name="grompp_nve",
                )
                self._mdrun("nve", "NVE")
                nve_summary = self._extract_nve_drift(warnings)

            update("محاكاة الإنتاج غير المقيدة (Parrinello-Rahman)...", 0.80)
            self._gmx(
                ["grompp", "-f", mdp["md.mdp"].name, "-c", "npt.gro", "-t", "npt.cpt",
                 "-p", staged.topol_top.name, "-n", "index.ndx", "-o", "prod_md.tpr",
                 "-maxwarn", str(cfg.grompp_maxwarn)],
                log_name="grompp_prod",
            )
            self._mdrun("prod_md", "Production")

            update("معالجة شروط الحدود الدورية (trjconv -pbc mol)...", 0.90)
            self._gmx(
                ["trjconv", "-s", "prod_md.tpr", "-f", "prod_md.xtc",
                 "-n", "index.ndx", "-o", "prod_md_center.xtc",
                 "-center", "-pbc", "mol", "-ur", "compact"],
                log_name="trjconv_center",
                input_text="Protein\nSystem\n",
            )
            self._gmx(
                ["trjconv", "-s", "prod_md.tpr", "-f", "prod_md_center.xtc",
                 "-n", "index.ndx", "-o", "prod_md_fit.xtc",
                 "-fit", "rot+trans"],
                log_name="trjconv_fit",
                input_text="Backbone\nProtein_Ligand\n",
            )
            self._gmx(
                ["trjconv", "-s", "prod_md.tpr", "-f", "prod_md_center.xtc",
                 "-n", "index.ndx", "-o", "trajectory_pbc_fit.pdb",
                 "-fit", "rot+trans"],
                log_name="trjconv_pdb",
                input_text="Backbone\nProtein_Ligand\n",
            )

            update("حساب RMSD وRMSF من المسار الفعلي...", 0.94)
            self._gmx(
                ["rms", "-s", "em.tpr", "-f", "prod_md_fit.xtc", "-n", "index.ndx",
                 "-o", "rmsd_ligand.xvg", "-tu", "ns"],
                log_name="rms_ligand",
                input_text="Backbone\nLigand_Heavy\n",
            )
            self._gmx(
                ["rmsf", "-s", "em.tpr", "-f", "prod_md_center.xtc", "-n", "index.ndx",
                 "-o", "rmsf_protein.xvg", "-res"],
                log_name="rmsf_protein",
                input_text="Protein\n",
            )

            update("تحليل الروابط الهيدروجينية (gmx hbond)...", 0.96)
            try:
                # GROMACS >= 2024 ships a new selection-based `hbond` that can
                # miss GAFF/AMBER ligand donors because it infers elements from
                # topology types.  The legacy group-based tool is exact and is
                # preferred when available; the new tool is a fallback.
                self._gmx(
                    ["hbond-legacy", "-f", "prod_md_center.xtc", "-s", "prod_md.tpr",
                     "-n", "index.ndx", "-num", "hbond_num.xvg"],
                    log_name="hbond",
                    input_text="Protein\nLigand\n",
                    timeout=3600,
                )
            except GromacsExecutionError:
                try:
                    self._gmx(
                        ["hbond", "-f", "prod_md_center.xtc", "-s", "prod_md.tpr",
                         "-n", "index.ndx", "-num", "hbond_num.xvg",
                         "-r", "Protein", "-t", "Ligand"],
                        log_name="hbond_new",
                        timeout=3600,
                    )
                except GromacsExecutionError as exc:
                    warnings.append(f"H-bond analysis failed: {str(exc)[:500]}")

            update("استخراج الطاقة والكثافة من ملف الطاقة (gmx energy)...", 0.98)
            energy_summary = self._extract_energy_summary(warnings)

            # Primary trajectory summary plus independent production replicates.
            replicate_stats.append(
                self._replicate_summary(
                    1, "prod_md", "rmsd_ligand.xvg", "hbond_num.xvg"
                )
            )
            for replicate_index in range(2, n_replicates + 1):
                update(
                    f"محاكاة إنتاج مستقلة (replicate {replicate_index}/{n_replicates})...",
                    0.98 + (0.01 * replicate_index / max(n_replicates, 1)),
                )
                replicate_stats.append(
                    self._run_extra_replicate(replicate_index, staged, mdp)
                )

        except Exception as exc:
            manifest = self._write_failed_manifest(str(exc), stage=failed_stage)
            raise SimulationError(str(exc)) from exc

        # Verify every claimed artifact exists.
        artifacts = {
            "complex_gro": staged.complex_gro.name,
            "complex_pdb": staged.complex_pdb.name,
            "topol_top": staged.topol_top.name,
            "ligand_itp": staged.ligand_itp.name,
            "ligand_posre": staged.ligand_posre.name,
            "index_ndx": "index.ndx",
            "solvated_gro": "solv_ions.gro",
            "em_gro": "em.gro",
            "production_tpr": "prod_md.tpr",
            "production_xtc": "prod_md.xtc",
            "trajectory_pbc": "prod_md_fit.xtc",
            "trajectory_pdb": "trajectory_pbc_fit.pdb",
            "rmsd_xvg": "rmsd_ligand.xvg",
            "rmsf_xvg": "rmsf_protein.xvg",
            "hbond_num_xvg": "hbond_num.xvg" if (out / "hbond_num.xvg").is_file() else None,
            "energy_edr": "prod_md.edr",
            "energy_xvg": "energy_production.xvg" if (out / "energy_production.xvg").is_file() else None,
        }
        missing = [name for name, fname in artifacts.items() if fname and not (out / fname).is_file()]
        if missing:
            message = (
                "GROMACS reported success but the following expected artifacts are missing: "
                + ", ".join(missing)
            )
            self._write_failed_manifest(message, stage="artifact_verification")
            raise SimulationError(message)

        update("اكتملت المحاكاة والتحليل بنجاح.", 1.00)

        def _aggregate(values: List[Optional[float]]) -> Dict[str, Any]:
            clean = [float(value) for value in values if value is not None]
            if not clean:
                return {"mean": None, "sd": None, "sem": None, "n": 0}
            mean = sum(clean) / len(clean)
            sd = (
                math.sqrt(sum((value - mean) ** 2 for value in clean) / (len(clean) - 1))
                if len(clean) > 1 else 0.0
            )
            return {
                "mean": mean,
                "sd": sd,
                "sem": sd / math.sqrt(len(clean)) if len(clean) > 1 else 0.0,
                "n": float(len(clean)),
            }

        replicate_statistics = {
            "n_replicates": len(replicate_stats),
            "rmsd_angstrom": _aggregate(
                [item["rmsd_mean_angstrom"] for item in replicate_stats]
            ),
            "hbond_count": _aggregate(
                [item["hbond_mean"] for item in replicate_stats]
            ),
        }

        manifest_data = {
            "run_id": cfg.run_id or out.name,
            "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "status": "completed",
            "engine_mode": "real_gromacs_execution",
            "software_versions": _software_versions(self.backend),
            "protocol_reference": "wiki/gromacs-protein-ligand-protocol.md",
            "engine": {
                "name": "GROMACS",
                "version": self.backend.version,
                "backend": self.backend.kind,
                "executable": self.backend.gmx_path,
            },
            "inputs_sha256": {
                name: _sha256_file(out / name)
                for name in (
                    staged.topol_top.name,
                    staged.complex_gro.name,
                    staged.ligand_itp.name,
                    "index.ndx",
                )
                if (out / name).is_file()
            },
            "topology": {
                "method": staged.method,
                "force_field": staged.force_field,
                "water_model": staged.water_model,
                "ligand_resname": staged.ligand_resname,
                "ligand_formal_charge": staged.ligand_charge,
                "protein_atoms": staged.protein_atom_count,
                "ligand_atoms": staged.ligand_atom_count,
                "water_validation": staged.water_validation,
                "retained_ions": staged.retained_ions,
                "retained_crystal_waters": staged.retained_waters,
            },
            "unit_cell": {
                "type": cfg.box_type,
                "distance_nm": cfg.box_distance_nm,
                "solute_bounding_diameter_nm": round(solute_diameter_nm, 3),
                "volume_savings_percent": savings_pct,
                "dodecahedron_volume_nm3": dodec_vol,
                "cubic_volume_nm3": cube_vol,
                "justification": (
                    f"Rhombic dodecahedral box is {savings_pct:.1f}% smaller than the equivalent "
                    "cubic box, reducing the number of explicit water molecules."
                ),
            },
            "parameters": {
                "temperature_k": cfg.temperature_k,
                "pressure_bar": cfg.pressure_bar,
                "timestep_ps": cfg.timestep_ps,
                "force_field": cfg.force_field,
                "water_model": cfg.water_model,
                "salt_molar": cfg.salt_concentration_molar,
                "nvt_duration_ps": cfg.nvt_duration_ps,
                "npt_duration_ps": cfg.npt_duration_ps,
                "production_duration_ns": cfg.prod_duration_ns,
                "random_seed": cfg.random_seed,
                "n_replicates": n_replicates,
                "analysis_burn_in_fraction": cfg.analysis_burn_in_fraction,
                "solvent_template": solvent_coordinate_file(cfg.water_model),
            },
            "index_groups": group_sizes,
            "analysis_reference": "em.tpr (energy-minimized starting structure)",
            "nve_check": nve_summary,
            "replicates": replicate_stats,
            "replicate_statistics": replicate_statistics,
            "replication_note": (
                "Replicate 1 is the primary trajectory. Extra replicates start from the "
                "same energy-minimized structure but use independent Maxwell-Boltzmann "
                "velocity seeds and independent NVT, NPT and production stages."
            ),
            "guardrails_enforced": {
                "ligand_heavy_restraints": True,
                "merged_thermostat_groups": "Protein_Ligand Water_and_ions",
                "berendsen_in_production": False,
                "production_barostat": "Parrinello-Rahman",
                "pbc_rewrapping_applied": True,
                "real_gromacs_execution": True,
            },
            "scripts": {k: str(v.name) for k, v in script_paths.items()},
            "artifacts": artifacts,
            "energy_summary": energy_summary,
            "limitations": [
                "Classical force field; no electronic polarisation or bond-breaking chemistry.",
                "A single production trajectory is not evidence of converged sampling.",
                f"Production length is {cfg.prod_duration_ns:g} ns; report only observables that are stable over that window.",
                "MM-PB(GB)SA is not computed by this pipeline; it requires gmx_MMPBSA on this trajectory.",
                "Ligand protonation/tautomer state is inherited from the docking preparation and is not re-optimised here.",
            ],
            "warnings": warnings,
        }

        manifest_path = out / "md_manifest.json"
        manifest_path.write_text(json.dumps(manifest_data, indent=2, ensure_ascii=False), encoding="utf-8")
        return manifest_data

    def _replicate_summary(
        self, index: int, prefix: str, rmsd_name: str, hbond_name: str
    ) -> Dict[str, Any]:
        rmsd_points = [
            (time_ns, value_nm * 10.0)
            for time_ns, value_nm in parse_xvg_file(self.output_dir / rmsd_name)
        ]
        hbond_points = parse_hbond_num_xvg(self.output_dir / hbond_name)
        if not rmsd_points:
            raise SimulationError(f"Replicate {index} has no RMSD data: {rmsd_name}")
        rmsd_summary = summarize_timeseries(
            rmsd_points, burn_in_fraction=self.config.analysis_burn_in_fraction
        )
        hbond_summary = summarize_timeseries(
            hbond_points, burn_in_fraction=self.config.analysis_burn_in_fraction
        )
        return {
            "replicate": index,
            "prefix": prefix,
            "rmsd_mean_angstrom": round(rmsd_summary["mean"], 4),
            "rmsd_sem_angstrom": round(rmsd_summary["sem"], 4),
            "rmsd_block_sem_angstrom": round(rmsd_summary["block_sem"], 4),
            "hbond_mean": round(hbond_summary["mean"], 4) if hbond_points else None,
            "hbond_sem": round(hbond_summary["sem"], 4) if hbond_points else None,
            "hbond_frames": int(hbond_summary["n"]),
            "frames": int(rmsd_summary["n"]),
        }

    def _run_extra_replicate(
        self, index: int, staged: TopologyResult, mdp: Dict[str, Path]
    ) -> Dict[str, Any]:
        """Run a fully independent NVT -> NPT -> production trajectory.

        Each replicate starts from the shared energy-minimized structure but gets
        a different Maxwell-Boltzmann seed for NVT, its own NPT equilibration,
        and its own production run.  This makes the replicate independent at the
        equilibration and production stages.
        """
        cfg = self.config
        seed = cfg.random_seed + index * 1000
        nvt_prefix = f"rep{index}_nvt"
        npt_prefix = f"rep{index}_npt"
        prefix = f"rep{index}_prod_md"

        nvt_mdp = self.output_dir / f"{nvt_prefix}.mdp"
        nvt_mdp.write_text(
            generate_mdp_nvt(
                duration_ps=cfg.nvt_duration_ps,
                temperature_k=cfg.temperature_k,
                timestep_ps=cfg.timestep_ps,
                force_field=cfg.force_field,
                random_seed=seed,
            ),
            encoding="utf-8",
            newline="\n",
        )
        self._gmx(
            ["grompp", "-f", nvt_mdp.name, "-c", "em.gro", "-r", "em.gro",
             "-p", staged.topol_top.name, "-n", "index.ndx",
             "-o", f"{nvt_prefix}.tpr", "-maxwarn", str(cfg.grompp_maxwarn)],
            log_name=f"grompp_{nvt_prefix}",
        )
        self._mdrun(nvt_prefix, f"Replicate {index} NVT")

        self._gmx(
            ["grompp", "-f", mdp["npt.mdp"].name, "-c", f"{nvt_prefix}.gro",
             "-t", f"{nvt_prefix}.cpt", "-r", f"{nvt_prefix}.gro",
             "-p", staged.topol_top.name, "-n", "index.ndx",
             "-o", f"{npt_prefix}.tpr", "-maxwarn", str(cfg.grompp_maxwarn)],
            log_name=f"grompp_{npt_prefix}",
        )
        self._mdrun(npt_prefix, f"Replicate {index} NPT")

        self._gmx(
            ["grompp", "-f", mdp["md.mdp"].name, "-c", f"{npt_prefix}.gro",
             "-t", f"{npt_prefix}.cpt", "-p", staged.topol_top.name,
             "-n", "index.ndx", "-o", f"{prefix}.tpr",
             "-maxwarn", str(cfg.grompp_maxwarn)],
            log_name=f"grompp_{prefix}",
        )
        self._mdrun(prefix, f"Replicate {index} production")
        self._gmx(
            ["trjconv", "-s", f"{prefix}.tpr", "-f", f"{prefix}.xtc",
             "-n", "index.ndx", "-o", f"{prefix}_center.xtc",
             "-center", "-pbc", "mol", "-ur", "compact"],
            log_name=f"trjconv_{prefix}_center",
            input_text="Protein\nSystem\n",
        )
        self._gmx(
            ["trjconv", "-s", f"{prefix}.tpr", "-f", f"{prefix}_center.xtc",
             "-n", "index.ndx", "-o", f"{prefix}_fit.xtc", "-fit", "rot+trans"],
            log_name=f"trjconv_{prefix}_fit",
            input_text="Backbone\nProtein_Ligand\n",
        )
        self._gmx(
            ["trjconv", "-s", f"{prefix}.tpr", "-f", f"{prefix}_center.xtc",
             "-n", "index.ndx", "-o", f"{prefix}_pbc_fit.pdb", "-fit", "rot+trans"],
            log_name=f"trjconv_{prefix}_pdb",
            input_text="Backbone\nProtein_Ligand\n",
        )
        self._gmx(
            ["rms", "-s", "em.tpr", "-f", f"{prefix}_fit.xtc", "-n", "index.ndx",
             "-o", f"{prefix}_rmsd_ligand.xvg", "-tu", "ns"],
            log_name=f"rms_{prefix}",
            input_text="Backbone\nLigand_Heavy\n",
        )
        try:
            self._gmx(
                ["hbond-legacy", "-f", f"{prefix}_center.xtc", "-s", f"{prefix}.tpr",
                 "-n", "index.ndx", "-num", f"{prefix}_hbond_num.xvg"],
                log_name=f"hbond_{prefix}",
                input_text="Protein\nLigand\n",
            )
        except GromacsExecutionError as exc:
            raise SimulationError(f"Replicate {index} H-bond analysis failed: {exc}") from exc
        return self._replicate_summary(
            index, prefix, f"{prefix}_rmsd_ligand.xvg", f"{prefix}_hbond_num.xvg"
        )

    def _extract_nve_drift(self, warnings: List[str]) -> Dict[str, Any]:
        """Extract total-energy drift from the short NVE segment."""
        try:
            self._gmx(
                ["energy", "-f", "nve.edr", "-o", "nve_energy.xvg"],
                log_name="energy_nve",
                input_text="Total-Energy\n\n",
                timeout=300,
            )
        except GromacsExecutionError as exc:
            warnings.append(f"NVE drift extraction failed: {str(exc)[:500]}")
            return {}

        # gmx energy was asked for Total-Energy only; retain its time column (ps).
        points = sorted(
            (time, value)
            for time, value in parse_xvg_file(self.output_dir / "nve_energy.xvg")
            if math.isfinite(time) and math.isfinite(value)
        )
        if len(points) < 2:
            warnings.append("NVE drift could not be computed: fewer than two total-energy frames.")
            return {}
        # Discard the first 10% (initial constraint/velocity transient) and fit a
        # linear slope rather than comparing two noisy endpoints.
        start = max(0, int(0.1 * len(points)))
        window = points[start:] if len(points) - start >= 2 else points
        n = len(window)
        mean_x = sum(time for time, _value in window) / n
        mean_y = sum(value for _time, value in window) / n
        denominator = sum((time - mean_x) ** 2 for time, _value in window)
        if denominator <= 0:
            warnings.append("NVE drift could not be computed: energy timestamps have no time span.")
            return {}
        slope = sum((time - mean_x) * (value - mean_y) for time, value in window) / denominator
        duration = points[-1][0] - points[0][0]
        try:
            n_atoms = len(parse_gro_atoms(self.output_dir / "solv_ions.gro")[0])
        except Exception:
            n_atoms = 0
        drift_per_atom = (slope / n_atoms) if n_atoms else None
        relative_percent_per_ps = abs(slope / mean_y * 100.0) if mean_y else None
        if not n_atoms:
            warnings.append("NVE drift per atom is unavailable: atom count could not be read.")
        summary = {
            "duration_ps": duration,
            "frames": len(points),
            "slope_kj_mol_per_ps": round(slope, 4),
            "delta_first_last_kj_mol": round(points[-1][1] - points[0][1], 3),
            "atoms": n_atoms,
            "drift_kj_mol_per_atom_per_ps": round(drift_per_atom, 6) if drift_per_atom is not None else None,
            "relative_drift_percent_per_ps": round(relative_percent_per_ps, 5) if relative_percent_per_ps is not None else None,
        }
        if ((drift_per_atom is not None and abs(drift_per_atom) > 0.01)
                or (relative_percent_per_ps is not None and relative_percent_per_ps > 0.1)):
            warnings.append(
                "NVE total-energy drift exceeds the screening threshold "
                f"(slope {slope:.5f} kJ/mol/ps); "
                "review timestep and nonbonded settings."
            )
        return summary

    def _extract_energy_summary(self, warnings: List[str]) -> Dict[str, Any]:
        out = self.output_dir
        selections = ["Temperature", "Pressure", "Density", "Potential", "Kinetic-En.", "Total-Energy"]
        expected_legends = {
            "Kinetic-En.": "Kinetic En.",
            "Total-Energy": "Total Energy",
        }
        try:
            self._gmx(
                ["energy", "-f", "prod_md.edr", "-o", "energy_production.xvg"],
                log_name="energy_production",
                input_text="\n".join(selections) + "\n\n",
                timeout=600,
            )
        except GromacsExecutionError as exc:
            warnings.append(f"Energy extraction failed: {str(exc)[:500]}")
            return {}

        columns = parse_gmx_energy_xvg(out / "energy_production.xvg")
        missing = [
            name for name in selections
            if expected_legends.get(name, name) not in columns
        ]
        if missing:
            warnings.append(
                "gmx energy did not return these requested terms: "
                + ", ".join(missing)
                + ". Only the available terms are reported."
            )
        summary: Dict[str, Any] = {}
        units = {
            "Temperature": "K",
            "Pressure": "bar",
            "Density": "kg/m^3",
            "Potential": "kJ/mol",
            "Kinetic En.": "kJ/mol",
            "Total Energy": "kJ/mol",
        }
        for name, values in columns.items():
            clean = [v for v in values if not math.isnan(v)]
            mean, std = _mean_std(clean)
            if mean is None:
                continue
            summary[name] = {"mean": mean, "std": std, "unit": units.get(name, ""), "n": len(clean)}
        return summary

    def _write_failed_manifest(self, error: str, stage: str) -> Dict[str, Any]:
        manifest = {
            "run_id": self.config.run_id or self.output_dir.name,
            "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "status": "failed",
            "engine_mode": "not_run",
            "failed_stage": stage,
            "error": error,
            "protocol_reference": "wiki/gromacs-protein-ligand-protocol.md",
            "engine": {
                "name": "GROMACS",
                "version": self.backend.version if self.backend else None,
                "backend": self.backend.kind if self.backend else None,
                "executable": self.backend.gmx_path if self.backend else None,
            },
        }
        try:
            (self.output_dir / "md_manifest.json").write_text(
                json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
            )
        except OSError:
            pass
        return manifest
