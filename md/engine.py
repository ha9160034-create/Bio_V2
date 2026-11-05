"""
GROMACS Molecular Dynamics Engine & Execution Layer
==================================================

Scientific Traceability & Wiki Documentation:
---------------------------------------------
- `wiki/gromacs-protein-ligand-protocol.md`: GROMACS Protein–Ligand Simulation Protocol (Justin Lemkul / LiveCoMS).
- `wiki/simulation-and-benchmark-best-practices.md`: LiveCoMS Foundations in Molecular Simulations.
- `wiki/open-force-field-and-ligand-parameterization.md`: Small-Molecule Force Field Parameterization.
- `wiki/force-fields.md`: Preserve the requested force-field version and parameter compatibility.
- `wiki/mm-pbsa-and-mm-gbsa.md`: Trajectory Requirements for MM-PB(GB)SA Rescoring.

Biophysical Guardrails (per AGENTS.md):
---------------------------------------
1. Rhombic dodecahedral unit cell default (`-bt dodecahedron -d 1.0`) to save ~23% solvent volume.
2. Restrain ligand heavy atoms during equilibration (k = 1000 kJ/(mol*nm^2)).
3. Merge Protein and Ligand into a single thermostat group (`tc-grps = Protein_Ligand Water_and_ions`).
4. NEVER use the Berendsen thermostat or barostat in production MD. Use Parrinello-Rahman barostat
   and Bussi velocity-rescaling (`v-rescale`) or Nosé-Hoover thermostat.
5. Always apply PBC re-wrapping (`trjconv -pbc mol -ur compact` and `-fit rot+trans`) before analysis.

Execution backends
------------------
GROMACS may be installed either:
  * natively (Linux/macOS/Windows), or
  * inside Windows Subsystem for Linux (WSL) on Windows.

This module detects both and normalises them behind :class:`GromacsBackend`.
No command is ever fabricated: if GROMACS cannot be found, callers receive
``None`` (discovery) or :class:`GromacsNotFoundError` (execution).
"""

from __future__ import annotations

import os
import sys
import re
import shlex
import shutil
import subprocess
import logging
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List, Sequence

logger = logging.getLogger("md.engine")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = PROJECT_ROOT / "tools"


class GromacsNotFoundError(RuntimeError):
    """Raised when a GROMACS command is requested but no installation exists."""


class GromacsExecutionError(RuntimeError):
    """Raised when a GROMACS command exits with a non-zero status."""


# ---------------------------------------------------------------------------
# Backend abstraction
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GromacsBackend:
    """A normalised handle to a working GROMACS installation."""

    kind: str                       # "native" | "wsl"
    executable: Tuple[str, ...]     # e.g. ("/usr/bin/gmx",) or ("wsl", "-e", "/usr/bin/gmx")
    gmx_path: str                   # path/command as reported to the user
    version: str = "unknown"
    data_prefix: str = ""

    @property
    def is_wsl(self) -> bool:
        return self.kind == "wsl"

    @property
    def label(self) -> str:
        where = "WSL" if self.is_wsl else "native"
        ver = self.version if self.version != "unknown" else "version unknown"
        return f"GROMACS {ver} ({where}: {self.gmx_path})"

    def command(self, args: Sequence[str]) -> List[str]:
        return [*self.executable, *[str(a) for a in args]]

    def to_backend_path(self, path: Path) -> str:
        """
        Translate a local filesystem path into the path syntax understood by
        the backend. On Windows + WSL this converts ``C:\\dir\\file`` to
        ``/mnt/c/dir/file``.
        """
        p = Path(path)
        try:
            p = p.resolve()
        except OSError:
            p = p.absolute()

        if not self.is_wsl:
            return str(p)

        drive = p.drive
        if not drive:
            # Already a POSIX path (or a relative path resolved on POSIX).
            return p.as_posix()

        letter = drive.rstrip(":").lower()
        tail = "/".join(p.parts[1:])
        return f"/mnt/{letter}/{tail}" if tail else f"/mnt/{letter}"

    def shell_command(self, script: str) -> List[str]:
        """Build a command that runs *script* through bash inside the backend."""
        if self.is_wsl:
            # executable[0] is wsl.exe; run bash -lc inside the distro.
            return [self.executable[0], "-e", "bash", "-lc", script]
        if sys.platform == "win32":
            raise GromacsExecutionError(
                "Shell commands (ACPYPE/OpenBabel) require a POSIX shell. "
                "Install WSL or run the platform on Linux/macOS."
            )
        return ["bash", "-lc", script]


_BACKEND_CACHE: Optional[GromacsBackend] = None
_BACKEND_PROBED = False
_FF_CACHE: Dict[str, List[str]] = {}


def _native_candidates() -> List[str]:
    names = ["gmx", "gmx_mpi"]
    if sys.platform == "win32":
        names = ["gmx.exe", "gmx_mpi.exe", "gmx", "gmx_mpi"]
    found: List[str] = []
    for name in names:
        which = shutil.which(name)
        if which:
            found.append(which)

    standard = [
        TOOLS_DIR / "gromacs" / "bin" / ("gmx.exe" if sys.platform == "win32" else "gmx"),
        TOOLS_DIR / "gmx" / ("gmx.exe" if sys.platform == "win32" else "gmx"),
        Path("C:/Program Files/GROMACS/bin/gmx.exe"),
        Path("C:/Program Files (x86)/GROMACS/bin/gmx.exe"),
        Path("/usr/local/gromacs/bin/gmx"),
        Path("/usr/bin/gmx"),
    ]
    for cand in standard:
        if cand.is_file():
            found.append(str(cand))
    return found


def _probe_native() -> Optional[GromacsBackend]:
    env_bin = os.environ.get("GMX_BIN") or os.environ.get("GROMACS_BIN")
    if env_bin:
        parts = shlex.split(env_bin, posix=(os.name != "nt"))
        if parts and Path(parts[0]).is_file():
            exe = str(Path(parts[0]).resolve())
            return GromacsBackend("native", (exe,), exe)
        if parts and parts[0].lower() in ("wsl", "wsl.exe"):
            wsl_exe = shutil.which(parts[0]) or parts[0]
            return GromacsBackend("wsl", (wsl_exe, *parts[1:]), " ".join(parts))
        if parts and shutil.which(parts[0]):
            return GromacsBackend("native", tuple(parts), parts[0])

    for cand in _native_candidates():
        return GromacsBackend("native", (cand,), cand)
    return None


def _probe_wsl(timeout: int = 25) -> Optional[GromacsBackend]:
    if sys.platform != "win32":
        return None
    wsl_exe = shutil.which("wsl") or shutil.which("wsl.exe")
    if not wsl_exe:
        return None
    try:
        res = subprocess.run(
            [wsl_exe, "-e", "bash", "-lc",
             "command -v gmx 2>/dev/null || command -v gmx_mpi 2>/dev/null"],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except Exception as exc:  # noqa: BLE001 - WSL may be absent/disabled
        logger.debug("WSL GROMACS probe failed: %s", exc)
        return None

    if res.returncode != 0:
        return None
    candidate = ""
    for line in res.stdout.splitlines():
        line = line.strip()
        if line.startswith("/") and " " not in line:
            candidate = line
            break
    if not candidate:
        return None
    return GromacsBackend("wsl", (wsl_exe, "-e", candidate), candidate)


def _enrich_backend(backend: GromacsBackend, timeout: int = 25) -> GromacsBackend:
    try:
        res = subprocess.run(
            backend.command(["--version"]),
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(PROJECT_ROOT),
        )
    except Exception as exc:  # noqa: BLE001
        logger.debug("Unable to query GROMACS version: %s", exc)
        return backend

    output = (res.stdout or "") + "\n" + (res.stderr or "")
    if res.returncode != 0:
        return backend

    version = "unknown"
    prefix = backend.data_prefix
    for raw in output.splitlines():
        line = raw.strip()
        if line.startswith("GROMACS version:"):
            version = line.split(":", 1)[1].strip()
        elif line.startswith("Data prefix:"):
            prefix = line.split(":", 1)[1].strip()
    return replace(backend, version=version, data_prefix=prefix)


def detect_gromacs(refresh: bool = False) -> Optional[GromacsBackend]:
    """
    Detect a working GROMACS installation.

    Order: explicit environment variable -> native PATH/standard prefixes ->
    WSL (Windows only). The result is cached because probing WSL costs ~1 s.
    Returns ``None`` when no usable installation exists.
    """
    global _BACKEND_CACHE, _BACKEND_PROBED
    if _BACKEND_PROBED and not refresh:
        return _BACKEND_CACHE

    backend = _probe_native()
    if backend is None:
        backend = _probe_wsl()
    if backend is not None:
        backend = _enrich_backend(backend)
        if backend.version == "unknown":
            logger.warning("GROMACS executable found but version query failed: %s", backend.gmx_path)
    _BACKEND_CACHE = backend
    _BACKEND_PROBED = True
    return backend


def find_gmx_path() -> Optional[Path]:
    """
    Backwards-compatible native-only discovery.

    Returns a local :class:`Path` for a native executable. WSL installations
    are reported through :func:`detect_gromacs` instead (their executable is
    ``wsl.exe``, not ``gmx``).
    """
    backend = _probe_native()
    if backend is None:
        return None
    return Path(backend.gmx_path)


def get_gmx_version(gmx_bin: Optional[Path] = None) -> Dict[str, Any]:
    """Retrieve GROMACS version and build capabilities (GPU, SIMD, MPI)."""
    if gmx_bin is not None:
        backend: Optional[GromacsBackend] = GromacsBackend(
            "native", (str(gmx_bin),), str(gmx_bin)
        )
        if backend is not None:
            backend = _enrich_backend(backend)
    else:
        backend = detect_gromacs()

    if backend is None:
        return {
            "available": False,
            "version": "Not installed",
            "path": None,
            "backend": None,
            "data_prefix": None,
            "gpu_acceleration": False,
            "simd": "None",
            "full_banner": "GROMACS binary was not detected natively or through WSL.",
        }

    try:
        res = subprocess.run(
            backend.command(["--version"]),
            capture_output=True,
            text=True,
            timeout=20,
            cwd=str(PROJECT_ROOT),
        )
        output = (res.stdout or "") + "\n" + (res.stderr or "")
        version_str = backend.version
        gpu_support = "GPU support:        disabled" not in output and (
            "CUDA" in output or "OpenCL" in output or "SYCL" in output
        )
        simd_info = "Default"
        for line in output.splitlines():
            line_str = line.strip()
            if line_str.startswith("GROMACS version:"):
                version_str = line_str.split(":", 1)[1].strip()
            elif "SIMD instructions:" in line_str:
                simd_info = line_str.split(":", 1)[1].strip()

        return {
            "available": True,
            "version": version_str,
            "path": backend.gmx_path,
            "backend": backend.kind,
            "data_prefix": backend.data_prefix,
            "gpu_acceleration": gpu_support,
            "simd": simd_info,
            "full_banner": output[:800],
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "available": False,
            "version": f"Error querying gmx ({exc})",
            "path": backend.gmx_path,
            "backend": backend.kind,
            "data_prefix": backend.data_prefix,
            "gpu_acceleration": False,
            "simd": "None",
            "full_banner": str(exc),
        }


# ---------------------------------------------------------------------------
# Force-field / water-model discovery
# ---------------------------------------------------------------------------

def _native_top_dir(backend: GromacsBackend) -> Optional[Path]:
    """Best-effort GROMACS ``share/gromacs/top`` directory for native installs."""
    if backend.is_wsl:
        return None
    prefix = Path(backend.data_prefix) if backend.data_prefix else None
    candidates: List[Path] = []
    if prefix:
        candidates.append(prefix / "share" / "gromacs" / "top")
        candidates.append(prefix / "gromacs" / "share" / "gromacs" / "top")
    exe = Path(backend.gmx_path)
    if exe.is_file():
        candidates.append(exe.parent.parent / "share" / "gromacs" / "top")
        candidates.append(exe.parent / "share" / "gromacs" / "top")
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return None


def list_force_fields(backend: Optional[GromacsBackend] = None) -> List[str]:
    """List force-field directories available to the backend (without .ff)."""
    backend = backend or detect_gromacs()
    if backend is None:
        return []
    if backend.data_prefix in _FF_CACHE:
        return _FF_CACHE[backend.data_prefix]

    # Native Windows has no POSIX shell.  Discover the share directory directly
    # instead of failing through shell_command().
    if sys.platform == "win32" and not backend.is_wsl:
        top_dir = _native_top_dir(backend)
        names = sorted(p.name[:-3] for p in top_dir.glob("*.ff")) if top_dir else []
        _FF_CACHE[backend.data_prefix] = names
        return names

    prefix = backend.data_prefix or ("/usr" if backend.is_wsl else "")
    top_dir = f"{prefix}/share/gromacs/top" if prefix else ""
    if not top_dir:
        return []

    script = f"ls -1d {shlex.quote(top_dir)}/*.ff 2>/dev/null | sed 's#.*/##; s/\\.ff$//' | sort"
    try:
        res = subprocess.run(
            backend.shell_command(script),
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(PROJECT_ROOT),
        )
        names = [ln.strip() for ln in res.stdout.splitlines() if ln.strip()]
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not list force fields: %s", exc)
        names = []

    _FF_CACHE[backend.data_prefix] = names
    return names


def resolve_force_field(backend: GromacsBackend, requested: str) -> str:
    """
    Map a human-readable force-field request to an installed GROMACS name.

    Raises:
        GromacsExecutionError: when the requested family is not installed.
    """
    available = list_force_fields(backend)
    key = (requested or "").strip().lower()

    aliases = {
        "amber": "amber99sb-ildn",
        "amber99sb-ildn": "amber99sb-ildn",
        "charmm36": "charmm36",
        "charmm36/cgenff": "charmm36",
        "charmm36-feb2026_cgenff-5.0": "charmm36",
        # OpenFF is not a GROMACS force-field directory and this deployment does
        # not implement a SMIRNOFF -> GROMACS export.  Do not silently map it to
        # Amber; let resolve_force_field raise an explicit error.
    }
    wanted = aliases.get(key, key)
    if wanted in available:
        return wanted

    if wanted == "charmm36":
        for cand in ("charmm36-jul2022", "charmm36-feb2022", "charmm36-jul2021", "charmm36"):
            if cand in available:
                return cand
        raise GromacsExecutionError(
            "CHARMM36 is not installed for this GROMACS backend. "
            f"Installed force fields: {', '.join(available) or 'none'}. "
            "Install charmm36.ff (and CGenFF parameters) or use Amber99SB-ILDN/GAFF2."
        )

    if wanted in ("amber99sb-ildn", "amber99sb", "amber14sb", "amber"):
        raise GromacsExecutionError(
            f"Requested force field '{requested}' is not installed for this GROMACS backend. "
            f"Installed force fields: {', '.join(available) or 'none'}."
        )

    raise GromacsExecutionError(
        f"Force field '{requested}' is not installed. "
        f"Installed force fields: {', '.join(available) or 'none'}."
    )


_WATER_ALIASES = {
    "tip3p": "tip3p",
    "tip3p (charmm/amber modified)": "tip3p",
    "spc": "spc",
    "spc/e": "spc",
    "spce": "spc",
    "tip4p": "tip4p",
    "tip4pew": "tip4pew",
    "tip4p-ew": "tip4pew",
    "tips3p": "tips3p",
}


def resolve_water_model(requested: str) -> str:
    return _WATER_ALIASES.get((requested or "").strip().lower(), (requested or "tip3p").strip().lower())


def solvent_coordinate_file(water_model: str) -> str:
    """Return the GROMACS ``-cs`` solvent template matching the water model.

    A 3-site template (``spc216.gro``) must not be used for a 4-site model such
    as TIP4P-Ew; mixing them creates a physically inconsistent solvent box.
    """
    model = resolve_water_model(water_model)
    if model in {"tip4p", "tip4pew", "tip5p"}:
        return "tip5p.gro" if model == "tip5p" else "tip4p.gro"
    return "spc216.gro"


_TOOL_CACHE: Dict[str, Dict[str, Any]] = {}


def check_backend_tools(backend: Optional[GromacsBackend] = None, refresh: bool = False) -> Dict[str, Any]:
    """
    Preflight the tools required by the topology route:
    ``gmx``, ``acpype``, ``antechamber`` and ``obabel``.

    Returns a dict with booleans and the list of installed force fields.
    """
    backend = backend or detect_gromacs()
    if backend is None:
        return {
            "gmx": False,
            "acpype": False,
            "antechamber": False,
            "obabel": False,
            "force_fields": [],
            "backend": None,
        }

    key = backend.gmx_path + "|" + backend.kind
    if key in _TOOL_CACHE and not refresh:
        return _TOOL_CACHE[key]

    script = (
        "for t in acpype obabel; do "
        "if command -v $t >/dev/null 2>&1; then echo \"$t=1\"; else echo \"$t=0\"; fi; done; "
        "if command -v antechamber >/dev/null 2>&1; then echo \"antechamber=1\"; "
        "elif command -v acpype >/dev/null 2>&1; then "
        "AP=$(python3 -c 'import acpype,os;print(os.path.dirname(acpype.__file__))' 2>/dev/null); "
        "if [ -n \"$AP\" ] && [ -x \"$AP/amber_linux/bin/antechamber\" ]; then echo \"antechamber=1\"; "
        "else echo \"antechamber=0\"; fi; else echo \"antechamber=0\"; fi"
    )
    tools = {"acpype": False, "antechamber": False, "obabel": False}
    if sys.platform == "win32" and not backend.is_wsl:
        # Native Windows: no POSIX shell.  Use PATH discovery and let the
        # native ACPYPE launcher be executed directly by the topology layer.
        tools["acpype"] = shutil.which("acpype") is not None
        tools["antechamber"] = shutil.which("antechamber") is not None or tools["acpype"]
        tools["obabel"] = shutil.which("obabel") is not None
    else:
        try:
            code, stdout, stderr = run_backend_shell(script, cwd=PROJECT_ROOT, timeout=60, backend=backend)
            for line in (stdout or "").splitlines():
                if "=" in line:
                    name, value = line.split("=", 1)
                    if name in tools:
                        tools[name] = value.strip() == "1"
        except Exception as exc:  # noqa: BLE001
            logger.warning("Tool preflight failed: %s", exc)

    result = {
        "gmx": True,
        "acpype": tools["acpype"],
        "antechamber": tools["antechamber"],
        "obabel": tools["obabel"],
        "force_fields": list_force_fields(backend),
        "backend": backend.kind,
        "version": backend.version,
    }
    _TOOL_CACHE[key] = result
    return result


# ---------------------------------------------------------------------------
# Command execution
# ---------------------------------------------------------------------------

def _looks_like_local_abs_path(arg: str) -> bool:
    if not arg:
        return False
    if re.match(r"^[A-Za-z]:[\\/]", arg):
        return True
    return arg.startswith("\\\\")


def run_gmx_command(
    command_args: List[str],
    cwd: Path,
    input_text: Optional[str] = None,
    log_name: Optional[str] = None,
    timeout: Optional[int] = None,
    backend: Optional[GromacsBackend] = None,
) -> Tuple[int, str, str]:
    """
    Execute a GROMACS command with full logging and output capture.

    Absolute Windows paths in ``command_args`` are translated to WSL paths
    when the backend is WSL. Relative filenames are left untouched because the
    working directory is carried over by ``wsl.exe``.
    """
    backend = backend or detect_gromacs()
    if backend is None:
        raise GromacsNotFoundError(
            "GROMACS ('gmx') was not found natively or through WSL. "
            "Install GROMACS or set GMX_BIN. Simulation cannot proceed."
        )

    cwd = Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)

    translated_args: List[str] = []
    for arg in command_args:
        arg_s = str(arg)
        if _looks_like_local_abs_path(arg_s):
            translated_args.append(backend.to_backend_path(Path(arg_s)))
        else:
            translated_args.append(arg_s)

    cmd = backend.command(translated_args)
    logger.info("Running: %s (cwd=%s)", " ".join(cmd), cwd)

    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd),
            input=input_text,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise GromacsNotFoundError(f"Unable to launch GROMACS backend {backend.label}: {exc}") from exc

    if log_name:
        log_file = cwd / f"{log_name}.log"
        with open(log_file, "w", encoding="utf-8") as f:
            f.write(f"COMMAND: {' '.join(cmd)}\n")
            f.write(f"CWD (host): {cwd}\n")
            f.write(f"CWD (backend): {backend.to_backend_path(cwd)}\n")
            f.write(f"EXIT CODE: {result.returncode}\n\n")
            f.write("--- STDOUT ---\n")
            f.write(result.stdout or "")
            f.write("\n--- STDERR ---\n")
            f.write(result.stderr or "")

    if result.returncode != 0:
        error_msg = (result.stderr or "").strip() or (result.stdout or "").strip()
        tool = command_args[0] if command_args else "gmx"
        logger.error("GROMACS %s failed (%s): %s", tool, backend.label, error_msg)
        raise GromacsExecutionError(
            f"GROMACS command '{tool}' failed with exit code {result.returncode}.\n"
            f"Backend: {backend.label}\n{error_msg}"
        )

    return result.returncode, result.stdout or "", result.stderr or ""


def run_backend_shell(
    script: str,
    cwd: Path,
    input_text: Optional[str] = None,
    log_name: Optional[str] = None,
    timeout: Optional[int] = None,
    backend: Optional[GromacsBackend] = None,
) -> Tuple[int, str, str]:
    """
    Run a bash command inside the backend (used for ACPYPE/OpenBabel).

    On Windows+WSL the working directory is translated and the command is
    executed inside the Linux filesystem view of the project.
    """
    backend = backend or detect_gromacs()
    if backend is None:
        raise GromacsNotFoundError("GROMACS backend is required to run ACPYPE/OpenBabel.")

    cwd = Path(cwd)
    cwd.mkdir(parents=True, exist_ok=True)
    translated_cwd = backend.to_backend_path(cwd)
    full_script = f"cd {shlex.quote(translated_cwd)} && set -o pipefail && {script}"
    cmd = backend.shell_command(full_script)

    try:
        result = subprocess.run(
            cmd,
            cwd=str(cwd),
            input=input_text,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise GromacsNotFoundError(f"Unable to launch shell backend: {exc}") from exc

    if log_name:
        log_file = cwd / f"{log_name}.log"
        with open(log_file, "w", encoding="utf-8") as f:
            f.write(f"SCRIPT: {script}\n")
            f.write(f"CWD (backend): {translated_cwd}\n")
            f.write(f"EXIT CODE: {result.returncode}\n\n")
            f.write("--- STDOUT ---\n")
            f.write(result.stdout or "")
            f.write("\n--- STDERR ---\n")
            f.write(result.stderr or "")

    return result.returncode, result.stdout or "", result.stderr or ""


# ---------------------------------------------------------------------------
# Standalone reproducibility scripts
# ---------------------------------------------------------------------------

def generate_standalone_scripts(
    output_dir: Path,
    complex_gro: str = "complex.gro",
    topol_top: str = "topol.top",
    index_ndx: str = "index.ndx",
    box_type: str = "dodecahedron",
    box_distance: float = 1.0,
    production_prefix: str = "prod_md",
    ligand_resname: str = "LIG",
    water_model: str = "tip3p",
    n_replicates: int = 1,
    salt_concentration_molar: float = 0.15,
) -> Dict[str, Path]:
    """
    Generate self-contained reproducible Bash (`run_md.sh`) and Windows
    (`run_md.bat`) scripts implementing the Lemkul/LiveCoMS protein-ligand
    protocol.

    The index file is regenerated from the actual solvated system by the
    bundled `make_index.py` helper (never hard-coded atom counts).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    solvent_gro = solvent_coordinate_file(water_model)
    replicate_indices = [str(index) for index in range(2, max(1, int(n_replicates)) + 1)]
    bash_replicates = ""
    bat_replicates = ""
    if replicate_indices:
        bash_replicates = f"""
# Additional fully independent replicates (NVT -> NPT -> production)
for i in {' '.join(replicate_indices)}; do
    echo ">>> Running independent replicate $i..."
    $GMX grompp -f nvt_rep$i.mdp -c em.gro -r em.gro -p {topol_top} -n {index_ndx} -o nvt_rep$i.tpr -maxwarn 2
    $GMX mdrun -v -deffnm nvt_rep$i
    $GMX grompp -f npt.mdp -c nvt_rep$i.gro -t nvt_rep$i.cpt -r nvt_rep$i.gro -p {topol_top} -n {index_ndx} -o npt_rep$i.tpr -maxwarn 2
    $GMX mdrun -v -deffnm npt_rep$i
    $GMX grompp -f md.mdp -c npt_rep$i.gro -t npt_rep$i.cpt -p {topol_top} -n {index_ndx} -o rep$i_prod_md.tpr -maxwarn 2
    $GMX mdrun -v -deffnm rep$i_prod_md
done
"""
        bat_replicates = f"""
rem Additional fully independent replicates (NVT -> NPT -> production)
for %%i in ({' '.join(replicate_indices)}) do (
    echo [Replicate %%i] NVT...
    %GMX% grompp -f nvt_rep%%i.mdp -c em.gro -r em.gro -p {topol_top} -n {index_ndx} -o nvt_rep%%i.tpr -maxwarn 2
    if errorlevel 1 exit /b %errorlevel%
    %GMX% mdrun -v -deffnm nvt_rep%%i
    %GMX% grompp -f npt.mdp -c nvt_rep%%i.gro -t nvt_rep%%i.cpt -r nvt_rep%%i.gro -p {topol_top} -n {index_ndx} -o npt_rep%%i.tpr -maxwarn 2
    if errorlevel 1 exit /b %errorlevel%
    %GMX% mdrun -v -deffnm npt_rep%%i
    %GMX% grompp -f md.mdp -c npt_rep%%i.gro -t npt_rep%%i.cpt -p {topol_top} -n {index_ndx} -o rep%%i_prod_md.tpr -maxwarn 2
    if errorlevel 1 exit /b %errorlevel%
    %GMX% mdrun -v -deffnm rep%%i_prod_md
)
"""

    bash_script = f"""#!/usr/bin/env bash
# ==============================================================================
# GROMACS All-Atom Protein-Ligand Molecular Dynamics Simulation Pipeline
# Protocol: Lemkul & LiveCoMS Foundations (see wiki/gromacs-protein-ligand-protocol.md)
# ==============================================================================
set -euo pipefail

echo "=========================================================="
echo " Starting In-Silico Molecular Dynamics Protocol"
echo "=========================================================="

GMX="${{GMX:-gmx}}"
PYTHON="${{PYTHON:-python3}}"

# 1. Unit Cell Definition (Rhombic Dodecahedron saves ~23% solvent volume)
echo ">>> [Step 1/10] Setting up {box_type} unit cell (d = {box_distance} nm)..."
$GMX editconf -f {complex_gro} -o newbox.gro -bt {box_type} -d {box_distance}

# 2. Explicit Solvation
echo ">>> [Step 2/10] Solvating system with explicit water..."
$GMX solvate -cp newbox.gro -cs {solvent_gro} -p {topol_top} -o solv.gro

# 3. Charge Neutralization & Configured Salt Concentration
echo ">>> [Step 3/10] Adding ions ({salt_concentration_molar:g} M NaCl)..."
$GMX grompp -f ions.mdp -c solv.gro -p {topol_top} -o ions.tpr -maxwarn 2
echo SOL | $GMX genion -s ions.tpr -o solv_ions.gro -p {topol_top} -pname NA -nname CL -neutral -conc {salt_concentration_molar:g}

# 4. Index groups built from the ACTUAL solvated system (no guessed counts)
echo ">>> [Step 4/10] Building index groups (Protein_Ligand / Water_and_ions)..."
$PYTHON make_index.py solv_ions.gro {index_ndx} {ligand_resname}

# 5. Energy Minimization (Steepest Descents, Fmax < 1000 kJ/(mol*nm))
echo ">>> [Step 5/10] Running energy minimization (em.mdp)..."
$GMX grompp -f em.mdp -c solv_ions.gro -p {topol_top} -o em.tpr -maxwarn 2
$GMX mdrun -v -deffnm em

# 6. NVT Thermalization (Bussi velocity-rescaling, merged thermostat groups)
echo ">>> [Step 6/10] Running NVT equilibration with heavy-atom restraints..."
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p {topol_top} -n {index_ndx} -o nvt.tpr -maxwarn 2
$GMX mdrun -v -deffnm nvt

# 7. NPT Isobaric Equilibration (Berendsen allowed only during equilibration)
echo ">>> [Step 7/10] Running NPT equilibration..."
$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p {topol_top} -n {index_ndx} -o npt.tpr -maxwarn 2
$GMX mdrun -v -deffnm npt

# 8. NVE Drift Check (no thermostat/barostat, from equilibrated NPT)
echo ">>> [Step 8/10] Running NVE drift check..."
$GMX grompp -f nve.mdp -c npt.gro -t npt.cpt -p {topol_top} -n {index_ndx} -o nve.tpr -maxwarn 2
$GMX mdrun -v -deffnm nve

# 9. Unrestrained Production MD (Parrinello-Rahman barostat, NO Berendsen)
echo ">>> [Step 9/10] Launching unrestrained production MD ({production_prefix})..."
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p {topol_top} -n {index_ndx} -o {production_prefix}.tpr -maxwarn 2
$GMX mdrun -v -deffnm {production_prefix}

# 10. Periodic Boundary Condition (PBC) Post-Processing
echo ">>> [Step 10/10] Re-wrapping trajectory (trjconv -pbc mol -ur compact)..."
printf "Protein\\nSystem\\n" | $GMX trjconv -s {production_prefix}.tpr -f {production_prefix}.xtc -n {index_ndx} -o {production_prefix}_center.xtc -center -pbc mol -ur compact
printf "Backbone\\nProtein_Ligand\\n" | $GMX trjconv -s {production_prefix}.tpr -f {production_prefix}_center.xtc -n {index_ndx} -o {production_prefix}_fit.xtc -fit rot+trans

{bash_replicates}
echo "=========================================================="
echo " Simulation and PBC post-processing completed successfully!"
echo "=========================================================="
"""

    bat_script = f"""@echo off
rem ==============================================================================
rem GROMACS All-Atom Protein-Ligand Molecular Dynamics Simulation Pipeline (Windows)
rem Protocol: Lemkul & LiveCoMS Foundations (see wiki/gromacs-protein-ligand-protocol.md)
rem ==============================================================================

echo ==========================================================
echo  Starting In-Silico Molecular Dynamics Protocol
echo ==========================================================

set GMX=gmx
set PYTHON=python

echo [Step 1/10] Setting up {box_type} unit cell...
%GMX% editconf -f {complex_gro} -o newbox.gro -bt {box_type} -d {box_distance}
if errorlevel 1 exit /b %errorlevel%

echo [Step 2/10] Solvating system with explicit water...
%GMX% solvate -cp newbox.gro -cs {solvent_gro} -p {topol_top} -o solv.gro
if errorlevel 1 exit /b %errorlevel%

echo [Step 3/10] Adding ions ({salt_concentration_molar:g} M NaCl)...
%GMX% grompp -f ions.mdp -c solv.gro -p {topol_top} -o ions.tpr -maxwarn 2
echo SOL | %GMX% genion -s ions.tpr -o solv_ions.gro -p {topol_top} -pname NA -nname CL -neutral -conc {salt_concentration_molar:g}
if errorlevel 1 exit /b %errorlevel%

echo [Step 4/10] Building index groups from the solvated system...
%PYTHON% make_index.py solv_ions.gro {index_ndx} {ligand_resname}
if errorlevel 1 exit /b %errorlevel%

echo [Step 5/10] Running energy minimization...
%GMX% grompp -f em.mdp -c solv_ions.gro -p {topol_top} -o em.tpr -maxwarn 2
%GMX% mdrun -v -deffnm em
if errorlevel 1 exit /b %errorlevel%

echo [Step 6/10] Running NVT equilibration...
%GMX% grompp -f nvt.mdp -c em.gro -r em.gro -p {topol_top} -n {index_ndx} -o nvt.tpr -maxwarn 2
%GMX% mdrun -v -deffnm nvt
if errorlevel 1 exit /b %errorlevel%

echo [Step 7/10] Running NPT equilibration...
%GMX% grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p {topol_top} -n {index_ndx} -o npt.tpr -maxwarn 2
%GMX% mdrun -v -deffnm npt
if errorlevel 1 exit /b %errorlevel%

echo [Step 8/10] Running NVE drift check...
%GMX% grompp -f nve.mdp -c npt.gro -t npt.cpt -p {topol_top} -n {index_ndx} -o nve.tpr -maxwarn 2
%GMX% mdrun -v -deffnm nve
if errorlevel 1 exit /b %errorlevel%

echo [Step 9/10] Running production MD (Parrinello-Rahman barostat)...
%GMX% grompp -f md.mdp -c npt.gro -t npt.cpt -p {topol_top} -n {index_ndx} -o {production_prefix}.tpr -maxwarn 2
%GMX% mdrun -v -deffnm {production_prefix}
if errorlevel 1 exit /b %errorlevel%

echo [Step 10/10] PBC Re-wrapping (trjconv -pbc mol -ur compact)...
(echo Protein & echo System) | %GMX% trjconv -s {production_prefix}.tpr -f {production_prefix}.xtc -n {index_ndx} -o {production_prefix}_center.xtc -center -pbc mol -ur compact
(echo Backbone & echo Protein_Ligand) | %GMX% trjconv -s {production_prefix}.tpr -f {production_prefix}_center.xtc -n {index_ndx} -o {production_prefix}_fit.xtc -fit rot+trans

{bat_replicates}
echo ==========================================================
echo  Simulation completed!
echo ==========================================================
"""

    sh_path = output_dir / "run_md.sh"
    bat_path = output_dir / "run_md.bat"

    with open(sh_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(bash_script)
    try:
        os.chmod(sh_path, 0o755)
    except Exception:
        pass

    with open(bat_path, "w", encoding="utf-8", newline="\r\n") as f:
        f.write(bat_script)

    return {"bash": sh_path, "batch": bat_path}
