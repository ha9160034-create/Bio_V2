# RCSB PDB coordinate downloads and mutation dataset fetching
import gzip
import json
import re
import requests
from pathlib import Path
import streamlit as st


def validate_mutation_database(data) -> tuple[bool, str]:
    """Validate the local mutation dataset before it is trusted.

    Schema: a JSON object mapping 4-character PDB IDs to 4-character PDB IDs.
    Anything else is rejected rather than partially used.
    """
    if not isinstance(data, dict):
        return False, "Mutation database must be a JSON object."
    pattern = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
    for key, value in data.items():
        if not isinstance(key, str) or not pattern.match(key.strip()):
            return False, f"Invalid PDB ID key in mutation database: {key!r}."
        if not isinstance(value, str) or not pattern.match(value.strip()):
            return False, f"Invalid PDB ID value for {key!r}: {value!r}."
    return True, f"{len(data)} entries validated."

@st.cache_data(ttl=3600)
def fetch_deposited_pdb(pdb_id: str):
    """Download the coordinates deposited in the PDB entry (the asymmetric unit)."""
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        return response.text
    except requests.RequestException as error:
        st.error(f"Could not download deposited coordinates for {pdb_id}: {error}")
        return None

@st.cache_data(ttl=3600)
def fetch_biological_assembly(pdb_id: str, assembly_id: int = 1):
    """Download the official RCSB biological assembly coordinate file.

    This file already contains the correct symmetry-generated oligomer. It must
    be used instead of guessing assemblies from chain names or sequences.
    """
    pdb_id = (pdb_id or "").strip().upper()
    if not pdb_id or pdb_id == "NONE":
        return None

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb{assembly_id}.gz"
    try:
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        return gzip.decompress(response.content).decode("utf-8")
    except (requests.RequestException, OSError, UnicodeDecodeError) as error:
        st.warning(
            f"Biological Assembly {assembly_id} is unavailable for {pdb_id}. "
            f"The deposited coordinates will be used unchanged. Details: {error}"
        )
        return None

@st.cache_data(ttl=60)
def load_mutation_db():
    """Load the *local, version-controlled* mutation database.

    The previous implementation silently fell back to an unpinned GitHub URL,
    so the scientific data could change without any trace in the project.
    Remote fetching has been removed: the database must be reviewed and
    committed to ``data/mutation.json``.  A malformed file is rejected rather
    than partially trusted.
    """
    base_dir = Path(__file__).resolve().parents[1]

    candidates = [base_dir / "data" / "mutation.json", base_dir / "mutation.json"]
    for candidate in candidates:
        if not candidate.exists():
            continue
        try:
            with open(candidate, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception as exc:
            st.error(f"تعذر قراءة قاعدة الطفرات المحلية {candidate.name}: {exc}")
            return {}
        valid, message = validate_mutation_database(data)
        if not valid:
            st.error(f"قاعدة الطفرات المحلية غير صالحة ({candidate.name}): {message}")
            return {}
        return data

    st.warning(
        "لم يتم العثور على قاعدة الطفرات المحلية `data/mutation.json`. "
        "لن يتم جلبها من الإنترنت تلقائياً لأن المصدر يجب أن يكون مثبتاً ومراجعاً."
    )
    return {}
