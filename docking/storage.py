# Storage and Manifest Management
import os
import json
import uuid
import hashlib
import datetime
from pathlib import Path

import sys
import shutil

PROJECT_ROOT = Path(__file__).resolve().parents[1]
VINA_WINDOWS = PROJECT_ROOT / "tools" / "vina" / "vina.exe"
VINA_LINUX = PROJECT_ROOT / "tools" / "vina" / "vina"
BASE_RUNS_DIR = PROJECT_ROOT / "docking_runs"

DISCLAIMER_TEXT = (
    "تنبيه علمي: هذه النتائج تمثل فرضية ارتباط حاسوبية (In-Silico Docking Hypothesis) "
    "تهدف لترتيب الوضعيات والمركبات ضمن هذا البروتوكول المحدد. "
    "لا تُعد هذه الدرجات قياساً معملياً لطاقة الارتباط الحرة أو دليلاً على الفاعلية السريرية أو الجرعة العلاجية."
)

# Pinned Vina binaries (PLAN.md §1.2). The Linux hash is enforced: a container
# that cannot present this exact binary must refuse to dock, because the apt
# package (1.2.7) would silently change every score.
VINA_SHA256 = {
    "linux": "fa0126a28a9ea9162d1b161dfa92bc76e632416db28ca246278ea4b2dc6860cb",
    "win32": "e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5",
}
VINA_EXPECTED_VERSION = {
    "linux": "AutoDock Vina v1.2.5",
    "win32": "AutoDock Vina v1.2.7",  # dev-only host binary
}


class VinaBinaryError(RuntimeError):
    """Raised when the Vina binary is not the pinned, verified artifact."""


def _vina_platform_key() -> str:
    return "win32" if sys.platform == "win32" else "linux"


def verify_vina_binary(path: Path) -> dict:
    """Return provenance for a Vina binary, enforcing the pin on Linux."""
    digest = compute_file_sha256(Path(path))
    key = _vina_platform_key()
    expected = VINA_SHA256[key]
    match = digest == expected
    if not match and key == "linux":
        raise VinaBinaryError(
            "Vina binary hash mismatch.\n"
            f"  path     : {path}\n"
            f"  sha256   : {digest}\n"
            f"  expected : {expected}\n"
            "Refusing to dock: results would not be comparable to the validated baseline."
        )
    return {
        "vina_path": str(path),
        "vina_sha256": digest,
        "vina_sha256_expected": expected,
        "vina_sha256_match": match,
        "vina_platform": key,
        "vina_expected_version": VINA_EXPECTED_VERSION[key],
    }


def get_vina_path() -> Path:
    if sys.platform == "win32":
        if VINA_WINDOWS.is_file():
            candidate = VINA_WINDOWS
        else:
            fallback = shutil.which("vina")
            candidate = Path(fallback) if fallback else VINA_WINDOWS
    else:
        if VINA_LINUX.is_file():
            try:
                os.chmod(VINA_LINUX, os.stat(VINA_LINUX).st_mode | 0o755)
            except Exception:
                pass
            candidate = VINA_LINUX
        else:
            fallback = shutil.which("vina")
            if fallback:
                # Fail closed: the apt-provided Vina (1.2.7) must never be used
                # as a silent substitute for the validated 1.2.5 binary.
                raise VinaBinaryError(
                    f"Pinned Vina binary is missing at {VINA_LINUX}; "
                    f"found {fallback} on PATH instead. Refusing to dock."
                )
            return VINA_LINUX
    if candidate.is_file():
        verify_vina_binary(candidate)
    return candidate

def create_run_directory(prefix="run") -> tuple[str, Path]:
    BASE_RUNS_DIR.mkdir(parents=True, exist_ok=True)
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S")
    uid = uuid.uuid4().hex[:6]
    run_id = f"{prefix}_{ts}_{uid}"
    run_dir = BASE_RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_id, run_dir

def write_initial_manifest(run_dir: Path, data: dict) -> Path:
    initial_path = run_dir / "manifest.initial.json"
    with open(initial_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    return initial_path

def finalize_manifest(run_dir: Path, data: dict) -> Path:
    manifest_path = run_dir / "manifest.json"
    temp_path = run_dir / "manifest.json.tmp"
    with open(temp_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    temp_path.replace(manifest_path)
    return manifest_path

def load_manifest(manifest_path: Path) -> dict:
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)

def compute_file_sha256(filepath: Path) -> str:
    """Compute SHA-256 checksum for manifest provenance tracking."""
    if not filepath.exists() or not filepath.is_file():
        return ""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()


def collect_software_versions() -> dict:
    """Record the exact software stack used to produce a docking manifest."""
    import sys

    versions = {"python": sys.version.split()[0]}
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
    try:
        from docking.engine import get_vina_version

        versions["autodock_vina"] = get_vina_version()
    except Exception:
        versions["autodock_vina"] = "unknown"
    try:
        vina_path = get_vina_path()
        if vina_path.is_file():
            versions.update(verify_vina_binary(vina_path))
    except Exception as exc:
        versions["vina_verification"] = f"failed: {exc}"
    return versions


def vina_version_warning(shown_version: str = "") -> str:
    """Return a warning when the runtime Vina version differs from the pin.

    The Linux container enforces the hash in get_vina_path(), so this is mainly
    for the Windows development path (1.2.7) where a different binary is allowed
    but results must not be compared with the validated 1.2.5 baseline.
    """
    key = _vina_platform_key()
    expected = VINA_EXPECTED_VERSION[key]
    if not shown_version:
        try:
            from docking.engine import get_vina_version

            shown_version = get_vina_version()
        except Exception:
            return ""
    if expected and expected not in shown_version:
        return (
            f"Vina version mismatch: running '{shown_version}', validated baseline is "
            f"'{expected}'. Scores are not directly comparable to published runs."
        )
    return ""

MANIFEST_STATUSES = {
    "validation passed",
    "validation passed (Top-10 recovery)",
    "validation failed",
    "validation unavailable",
    "protocol validated on WT reference",
    "matched comparison—mutant arm (WT validation not transferable)",
    "not validated for affinity comparison",
    "exploratory—no reference validation",
    "exploratory—insufficient independent seeds",
    "exploratory—cross-protein docking",
    "run failed",
}

LEGACY_STATUS_ALIASES = {
    "exploratory": "exploratory—no reference validation",
}


STATUS_MAP = {s.lower(): s for s in MANIFEST_STATUSES}
STATUS_MAP.update({k.lower(): v for k, v in LEGACY_STATUS_ALIASES.items()})


def normalize_manifest_status(status: str) -> str:
    key = str(status).strip().lower()
    if key in STATUS_MAP:
        return STATUS_MAP[key]
    raise ValueError(f"Unsupported docking manifest status: {status!r}")


def build_and_save_manifest(
    run_dir: Path,
    run_id: str,
    status: str,
    target_meta: dict,
    binding_site_meta: dict,
    ligand_meta: dict,
    prep_meta: dict,
    engine_meta: dict,
    grid_meta: dict,
    results_meta: dict,
    warnings: list[str] = None,
    input_fingerprint: str = None
) -> Path:
    """Build and write the final manifest for a docking run."""
    run_dir = Path(run_dir)
    status = normalize_manifest_status(status)

    created_at_utc = target_meta.get("created_at_utc") or ""

    inputs_dict = {
        "sha256": {
            p.name: compute_file_sha256(p)
            for p in sorted(run_dir.iterdir())
            if p.is_file() and not p.name.startswith("manifest")
        }
    }
    if input_fingerprint:
        inputs_dict["fingerprint"] = input_fingerprint

    manifest_data = {
        "schema_version": 1,
        "run_id": run_id,
        "status": status,
        "created_at_utc": created_at_utc,
        "input_fingerprint": input_fingerprint,
        "target": target_meta,
        "binding_site": binding_site_meta,
        "ligand": ligand_meta,
        "receptor_preparation": prep_meta,
        "engine": engine_meta,
        "software_versions": collect_software_versions(),
        "grid": grid_meta,
        "inputs": inputs_dict,
        "results": results_meta,
        "warnings": warnings or [],
        "scientific_disclaimer": DISCLAIMER_TEXT
    }
    return finalize_manifest(run_dir, manifest_data)
