"""
RCSB PDB and Biological Assembly coordinate downloading service.
Independent module to avoid circular imports between b.py and docking modules.
"""
import gzip
import statistics
import xml.etree.ElementTree as ET

import requests

_AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",
}


def fetch_validation_xml(pdb_id: str, timeout: int = 20) -> bytes | None:
    """Download the official wwPDB validation XML (gzip) for an entry."""
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None
    url = f"https://files.rcsb.org/validation/download/{pdb_id.lower()}_validation.xml.gz"
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code != 200:
            return None
        return gzip.decompress(response.content)
    except Exception:
        return None


def parse_ligand_rscc(xml_bytes: bytes, ligand_code: str) -> dict:
    """Parse ligand RSCC values from a wwPDB validation XML report.

    RSCC is a model-to-density fit metric; it is not derived from docking and
    must not be presented as a docking score.  Multiple chains/instances are
    returned separately.
    """
    code = (ligand_code or "").strip().upper()
    if not code:
        return {"ligand_code": "", "ligand_rscc": None, "entries": []}
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        return {"ligand_code": code, "ligand_rscc": None, "entries": []}

    entries = []
    for element in root.iter("ModelledSubgroup"):
        if (element.get("resname") or "").strip().upper() != code:
            continue
        raw_rscc = element.get("rscc")
        try:
            rscc = float(raw_rscc) if raw_rscc not in (None, "", ".", "?") else None
        except ValueError:
            rscc = None
        entries.append({
            "chain": (element.get("chain") or "").strip(),
            "resnum": (element.get("resnum") or "").strip(),
            "icode": (element.get("icode") or "").strip(),
            "rscc": rscc,
            "occupancy": element.get("avgoccu"),
        })
    values = [entry["rscc"] for entry in entries if entry["rscc"] is not None]
    return {
        "ligand_code": code,
        "ligand_rscc": min(values) if values else None,
        "entries": entries,
    }


def parse_residue_rscc(xml_bytes: bytes) -> dict:
    """Parse protein-residue RSCC values from a wwPDB validation report.

    Multiple altlocs are collapsed conservatively to the lowest RSCC for each
    ``(chain, resnum, icode)`` key.  Waters, ligands and other HETATM groups are
    excluded from the protein summary.
    """
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        return {}
    residues: dict = {}
    for element in root.iter("ModelledSubgroup"):
        resname = (element.get("resname") or "").strip().upper()
        if resname not in _AA3:
            continue
        raw_rscc = element.get("rscc")
        try:
            rscc = float(raw_rscc) if raw_rscc not in (None, "", ".", "?") else None
        except ValueError:
            rscc = None
        if rscc is None:
            continue
        key = (
            (element.get("chain") or "").strip(),
            (element.get("resnum") or "").strip(),
            (element.get("icode") or "").strip(),
        )
        if key not in residues or rscc < residues[key]:
            residues[key] = rscc
    return residues


def summarize_residue_rscc(residues: dict) -> dict:
    """Conservative summary of protein-residue RSCC values."""
    values = list(residues.values())
    if not values:
        return {
            "n_residues": 0,
            "median_rscc": None,
            "min_rscc": None,
            "max_rscc": None,
            "n_below_0_8": 0,
            "fraction_below_0_8": None,
        }
    return {
        "n_residues": len(values),
        "median_rscc": round(float(statistics.median(values)), 4),
        "min_rscc": round(float(min(values)), 4),
        "max_rscc": round(float(max(values)), 4),
        "n_below_0_8": sum(1 for value in values if value < 0.8),
        "fraction_below_0_8": round(
            sum(1 for value in values if value < 0.8) / len(values), 4
        ),
    }


def fetch_ligand_rscc(pdb_id: str, ligand_code: str, timeout: int = 20) -> dict:
    """Download and parse the official wwPDB validation XML for one ligand."""
    code = (ligand_code or "").strip().upper()
    if not code:
        return {"ligand_code": "", "ligand_rscc": None, "entries": []}
    xml_bytes = fetch_validation_xml(pdb_id, timeout=timeout)
    if xml_bytes is None:
        return {"ligand_code": code, "ligand_rscc": None, "entries": []}
    return parse_ligand_rscc(xml_bytes, code)


def fetch_deposited_pdb_text(pdb_id: str, timeout: int = 60):
    """Download deposited coordinates (asymmetric unit) from RCSB."""
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.text
    except Exception:
        return None

def fetch_biological_assembly_text(pdb_id: str, assembly_id: int = 1, timeout: int = 60):
    """Download official RCSB biological assembly (.pdb1.gz)."""
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb{assembly_id}.gz"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return gzip.decompress(response.content).decode("utf-8")
    except Exception:
        return None


def fetch_entry_quality(
    pdb_id: str,
    ligand_code: str = "",
    timeout: int = 20,
    include_residue_rscc: bool = False,
) -> dict:
    """Fetch crystallographic quality metadata (resolution, R_free, method, RSCC).

    Returns an empty dict when metadata is unavailable.  Callers should treat
    missing values as "unknown", not as passing a quality threshold.
    """
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return {}
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code != 200:
            return {}
        data = response.json()
    except Exception:
        return {}

    entry_info = data.get("rcsb_entry_info") or {}
    refine = (data.get("refine") or [{}])[0]
    resolutions = entry_info.get("resolution_combined") or []
    resolution = float(resolutions[0]) if resolutions else None
    r_free = refine.get("ls_R_factor_R_free")
    r_work = refine.get("ls_R_factor_R_work")
    result = {
        "pdb_id": pdb_id,
        "title": (data.get("struct") or {}).get("title", ""),
        "experimental_method": entry_info.get("experimental_method", ""),
        "resolution_angstrom": resolution,
        "r_free": float(r_free) if r_free is not None else None,
        "r_work": float(r_work) if r_work is not None else None,
        "meets_resolution_target": resolution is not None and resolution <= 2.5,
        "meets_rfree_target": r_free is not None and float(r_free) < 0.25,
    }
    if ligand_code or include_residue_rscc:
        xml_bytes = fetch_validation_xml(pdb_id, timeout=timeout)
        if xml_bytes is not None:
            if ligand_code:
                result.update(parse_ligand_rscc(xml_bytes, ligand_code))
            if include_residue_rscc:
                residues = parse_residue_rscc(xml_bytes)
                result["residue_rscc_summary"] = summarize_residue_rscc(residues)
                result["residue_rscc_low"] = [
                    {
                        "chain": chain,
                        "resnum": resnum,
                        "icode": icode,
                        "rscc": round(rscc, 3),
                    }
                    for (chain, resnum, icode), rscc in sorted(
                        residues.items(), key=lambda item: item[1]
                    )
                    if rscc < 0.8
                ][:25]
    return result


def fetch_chem_comp_smiles(comp_id: str, timeout: int = 15) -> str | None:
    """Download verified reference SMILES for a Chemical Component ID from RCSB."""
    comp_id = (comp_id or "").strip().upper()
    if not comp_id:
        return None
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code == 200:
            data = response.json()
            descriptor = data.get("rcsb_chem_comp_descriptor") or {}
            smiles = descriptor.get("SMILES_stereo") or descriptor.get("SMILES")
            if smiles:
                return str(smiles).strip()
        return None
    except Exception:
        return None
