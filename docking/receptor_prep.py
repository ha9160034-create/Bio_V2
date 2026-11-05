# Receptor preparation for docking
import os
from pathlib import Path
from io import StringIO
from Bio.PDB import PDBParser, PDBIO
from docking.openbabel_io import read_file, write_file

_HET_ELEMENT_CODES = {"CL", "BR", "NA", "MG", "ZN", "FE", "CA", "MN", "CU", "K", "LI", "RB", "CS", "IOD"}
_HEME_RESNAMES = {"HEM", "HEC", "HEB", "HEA"}


def _infer_element_from_context(record: str, atom_name: str, resname: str) -> str:
    """Infer an element symbol without the classic ``CA`` = calcium trap.

    In protein ``ATOM`` records ``CA`` is the alpha carbon, never calcium.
    Two-letter elements are only accepted for ``HETATM`` records whose residue
    name is the matching ion (e.g. ``CA``/``CL``/``NA``) or for the iron of a
    heme group. Everything else falls back to the first letter of the atom name.
    """
    name = (atom_name or "").strip()
    alpha = "".join(ch for ch in name if ch.isalpha())
    if not alpha:
        return "C"
    upper = alpha.upper()
    res = (resname or "").strip().upper()

    if record == "HETATM":
        if res in _HET_ELEMENT_CODES and upper == res:
            return res
        if res in _HEME_RESNAMES and upper in {"FE", "NA", "MG", "CA"}:
            return upper
    elif res == "MSE" and upper == "SE":
        return "SE"

    return alpha[0].upper()


def _format_pdb_element_columns(pdb_path: Path) -> None:
    """
    Standardize ATOM and HETATM records to strictly adhere to PDB format specifications:
    columns 77-78 (1-indexed) must contain the element symbol right-justified.
    BioPython PDBIO writes the element shifted by 1 column (columns 76-77, left-justified),
    which causes OpenBabel to fail element parsing and misidentify atom names like CMA as Curium (Cm),
    or splits two-letter elements like FE into 'F E'.
    """
    formatted_lines = []
    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                line_str = line.rstrip("\r\n")
                atom_name = line_str[12:16].strip() if len(line_str) >= 16 else ""
                resname = line_str[17:20].strip() if len(line_str) >= 20 else ""
                record = line_str[:6].strip()
                # Inspect columns 73-80 to extract the raw element written by BioPython
                tail = line_str[72:80] if len(line_str) >= 80 else (line_str[72:] if len(line_str) > 72 else "")
                raw_elem = "".join(ch for ch in tail if ch.isalpha()).strip()
                if not raw_elem or len(raw_elem) > 2:
                    raw_elem = _infer_element_from_context(record, atom_name, resname)
                elif len(raw_elem) == 2 and raw_elem.upper() not in _HET_ELEMENT_CODES and record == "ATOM":
                    raw_elem = _infer_element_from_context(record, atom_name, resname)
                # Columns 1-76 must be exactly 76 characters long without stray letters at col 76 (char index 75)
                prefix = line_str[:72].ljust(76)
                elem_rjust = raw_elem.upper().rjust(2) if len(raw_elem) <= 2 else raw_elem[:2].upper().rjust(2)
                formatted_lines.append(f"{prefix}{elem_rjust}\n")
            else:
                formatted_lines.append(line)
    with open(pdb_path, "w", encoding="utf-8") as f:
        f.writelines(formatted_lines)


def _sanitize_pdb_for_docking(pdb_string: str) -> tuple[str, dict]:
    """Keep the first model and the highest-occupancy alternate location.

    Vina/OpenBabel do not consistently resolve multi-model or alternate-location
    records.  Passing duplicate atoms into a rigid receptor can create spurious
    contacts or ambiguous atom typing, so the choice is made explicitly and
    recorded in the preparation report.
    """
    lines = pdb_string.splitlines()
    model_lines: list[str] = []
    models_seen = 0
    in_first_model = True
    for line in lines:
        record = line[:6].strip()
        if record == "MODEL":
            models_seen += 1
            if models_seen > 1:
                in_first_model = False
            continue
        if record == "ENDMDL":
            if models_seen <= 1:
                in_first_model = False
            continue
        if models_seen > 1 and not in_first_model:
            continue
        model_lines.append(line)

    grouped: dict[tuple, tuple[float, int, str]] = {}
    order: list[tuple] = []
    removed_altloc = 0
    output: list[str] = []
    for line in model_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            output.append(line)
            continue
        altloc = line[16:17]
        if altloc in (" ", ""):
            output.append(line)
            continue
        removed_altloc += 1
        key = (
            line[:6], line[17:20], line[21:22], line[22:27], line[12:16]
        )
        try:
            occupancy = float(line[54:60])
        except ValueError:
            occupancy = 0.0
        previous = grouped.get(key)
        score = (occupancy, 1 if altloc in ("A", "1") else 0)
        if previous is None:
            grouped[key] = (occupancy, score[1], line)
            order.append(key)
        elif score > (previous[0], previous[1]):
            grouped[key] = (occupancy, score[1], line)

    # Re-insert chosen altloc records at their first occurrence position.
    chosen_lines = {key: value[2] for key, value in grouped.items()}
    emitted: set[tuple] = set()
    sanitized: list[str] = []
    for line in model_lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            sanitized.append(line)
            continue
        key = (
            line[:6], line[17:20], line[21:22], line[22:27], line[12:16]
        )
        if line[16:17] in (" ", ""):
            if key not in emitted:
                sanitized.append(line)
                emitted.add(key)
            continue
        if key in emitted:
            continue
        chosen = chosen_lines.get(key, line)
        if len(chosen) >= 17 and chosen[16:17] not in (" ", ""):
            chosen = chosen[:16] + " " + chosen[17:]
        sanitized.append(chosen)
        emitted.add(key)

    stats = {
        "models_dropped": max(0, models_seen - 1),
        "altloc_records_resolved": removed_altloc,
        "model_policy": "first_model_only" if models_seen > 1 else "single_model",
        "altloc_policy": "highest_occupancy_then_A" if removed_altloc else "none_present",
    }
    return "\n".join(sanitized) + "\n", stats


def prepare_receptor(
    pdb_string: str,
    output_dir: Path,
    keep_waters: bool = False,
    reference_ligand_resname: str = None
) -> tuple[Path, Path, dict]:
    output_dir = Path(output_dir)
    pdb_string, sanitizer_stats = _sanitize_pdb_for_docking(pdb_string)
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("target", StringIO(pdb_string))

    removed_waters = 0
    removed_ligands = []
    retained_cofactors = []
    reference_ligand_removed = False
    ref_ligand_lines = []

    ref_upper = (reference_ligand_resname or "").strip().upper()

    # 1. Preserve the reference ligand in reference_ligand.pdb before receptor removal
    if ref_upper:
        for line in pdb_string.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                rname = line[17:20].strip().upper()
                if rname == ref_upper:
                    ref_ligand_lines.append(line)
        if ref_ligand_lines:
            ref_lig_path = output_dir / "reference_ligand.pdb"
            with open(ref_lig_path, "w", encoding="utf-8") as f:
                f.write("\n".join(ref_ligand_lines) + "\n")

    # 2. Separate protein, ligands, waters, and cofactors with full provenance (resname, chain, resnum)
    for model in struct:
        for chain in model:
            to_detach = []
            for res in chain:
                res_id = res.id[0]
                res_name = res.resname.strip().upper()
                res_num = int(res.id[1]) if isinstance(res.id[1], int) else str(res.id[1])
                chain_id = str(chain.id)
                res_entry = {"resname": res_name, "chain": chain_id, "resnum": res_num}

                if res_id.startswith("W"):
                    if not keep_waters:
                        to_detach.append(res.id)
                        removed_waters += 1
                elif res_id != " ":
                    if ref_upper and res_name == ref_upper:
                        to_detach.append(res.id)
                        reference_ligand_removed = True
                    elif res_name in ["HEM", "MG", "ZN", "CA", "FE", "MN"]:
                        retained_cofactors.append(res_entry)
                    else:
                        to_detach.append(res.id)
                        removed_ligands.append(res_entry)
            for rid in to_detach:
                chain.detach_child(rid)

    # 3. Save immutable receptor clean PDB
    clean_pdb_path = output_dir / "receptor_clean.pdb"
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(clean_pdb_path.resolve()))
    _format_pdb_element_columns(clean_pdb_path)

    # 4. Convert to standard rigid PDBQT
    clean_pdbqt_path = output_dir / "receptor.pdbqt"
    mol = read_file(clean_pdb_path, "pdb")
    # OpenBabel adds hydrogens with pH-aware protonation (default pH 7.4). The
    # PDBQT writer merges nonpolar hydrogens into their heavy atoms, so the
    # written receptor contains AutoDock polar hydrogens (HD) only.
    mol.AddHydrogens(False, True, 7.4)
    write_file(mol, clean_pdbqt_path, "pdbqt", rigid_receptor=True)

    # 5. Sanitize rigid receptor PDBQT:
    # - Strictly remove ROOT/BRANCH/TORSDOF tags
    # - Fix any unsupported or misformatted AutoDock atom types (e.g. Cm from heme CMA)
    # - Guarantee that every ATOM/HETATM line has valid non-whitespace columns 78-79 (len >= 79)
    #   to prevent AutoDock Vina parse_pdbqt.cpp(69) assertion crash.
    lines = []
    with open(clean_pdbqt_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
                continue
            if line.startswith(("ATOM  ", "HETATM")):
                clean_line = line.rstrip("\r\n")
                if len(clean_line) < 77:
                    clean_line = clean_line.ljust(77)
                ad_type = clean_line[77:79].strip() if len(clean_line) >= 79 else clean_line[77:].strip()
                if ad_type == "Cm":
                    ad_type = "C"
                if not ad_type:
                    atom_name = clean_line[12:16].strip() if len(clean_line) >= 16 else ""
                    res_name = clean_line[17:20].strip() if len(clean_line) >= 20 else ""
                    inferred_el = _infer_element_from_context(clean_line[:6].strip(), atom_name, res_name)
                    if inferred_el in ("FE", "Fe"):
                        ad_type = "Fe"
                    elif inferred_el in ("ZN", "Zn"):
                        ad_type = "Zn"
                    elif inferred_el in ("MG", "Mg"):
                        ad_type = "Mg"
                    elif inferred_el in ("CA", "Ca"):
                        ad_type = "Ca"
                    elif inferred_el in ("MN", "Mn"):
                        ad_type = "Mn"
                    elif inferred_el in ("CL", "Cl"):
                        ad_type = "Cl"
                    elif inferred_el in ("BR", "Br"):
                        ad_type = "Br"
                    elif inferred_el in ("C", "A", "N", "NA", "O", "OA", "S", "SA", "P", "F", "I", "HD"):
                        ad_type = inferred_el
                    else:
                        ad_type = inferred_el[0].upper()
                type_formatted = ad_type.ljust(2) if len(ad_type) == 1 else ad_type[:2]
                clean_line = clean_line[:77].ljust(77) + type_formatted + "\n"
                lines.append(clean_line)
            else:
                lines.append(line)
    with open(clean_pdbqt_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    prep_report = {
        "clean_pdb_path": str(clean_pdb_path.name),
        "clean_pdbqt_path": str(clean_pdbqt_path.name),
        "reference_ligand_saved": bool(ref_ligand_lines),
        "reference_ligand_removed": reference_ligand_removed,
        "waters_policy": "Bulk waters removed (standard rigid-receptor docking model assumption)" if not keep_waters else "Retained waters",
        "waters_removed_count": removed_waters,
        "ligands_removed": removed_ligands,
        "cofactors_retained": retained_cofactors,
        "cofactor_parameterization_note": "Retained cofactors are processed with standard Gasteiger charges; no advanced force-field parameterization is claimed.",
        "hydrogen_policy": (
            "All hydrogens were added with OpenBabel pH-aware protonation (pH 7.4); "
            "the PDBQT writer merges nonpolar hydrogens, leaving AutoDock polar HD atoms."
        ),
        "all_hydrogens_added_with_ph_correction": True,
        "protonation_ph": 7.4,
        "pdbqt_atom_typing": "AutoDock atom types via OpenBabel (nonpolar H merged)",
        **sanitizer_stats,
    }
    return clean_pdb_path, clean_pdbqt_path, prep_report
