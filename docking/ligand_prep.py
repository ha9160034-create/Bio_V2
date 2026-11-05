# Ligand preparation using RDKit and OpenBabel
"""Small-molecule preparation for AutoDock Vina.

Scientific provenance
---------------------
* Protonation is applied explicitly at the requested pH with OpenBabel's
  pH-aware hydrogen model (``AddHydrogens(correctForPH=True, pH=...)``), then a
  fresh 3D conformer is generated with RDKit and AM1-BCC/GAFF2 parameters are
  produced later by ACPYPE.  The reported ``protonation_applied`` flag is only
  ``True`` when this step actually ran.
* Tautomers are **not** enumerated.  The chosen protomer is reported verbatim so
  that a reviewer can see exactly which chemical species was docked.
* A failure to apply pH protonation is fail-closed: the pipeline refuses to
  claim a pH-7.4 preparation when it could not perform one.
"""
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from docking.openbabel_io import read_file, write_file, write_text


def _build_3d_conformer(mol: Chem.Mol) -> Chem.Mol:
    """Add hydrogens and generate an MMFF-optimised 3D conformer."""
    mol = Chem.AddHs(mol)
    embed_code = AllChem.EmbedMolecule(mol, randomSeed=42)
    if embed_code != 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    return mol


def _apply_ph_protonation(smiles: str, ph: float) -> tuple[str, str]:
    """Return ``(protomer_smiles, method)`` at the requested pH.

    OpenBabel's pH model is a rule-based pKa correction, not a full tautomer
    enumerator.  It is used because it is the only pH-aware protonation engine
    available in this deployment; the exact protomer is recorded in the report.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: '{smiles}'")

    # A hydrogen-suppressed molecule gives OpenBabel full control over which
    # protons are present at the requested pH.
    heavy = Chem.RemoveHs(mol)
    heavy_block = Chem.MolToMolBlock(heavy)
    from docking.openbabel_io import read_text

    ob_mol = read_text(heavy_block, "mol")
    ob_mol.DeleteHydrogens()
    if not ob_mol.AddHydrogens(False, True, float(ph)):
        raise ValueError("OpenBabel could not add pH-corrected hydrogens.")
    protomer = write_text(ob_mol, "can").splitlines()[0].split()[0].strip()
    if not protomer:
        raise ValueError("OpenBabel returned an empty canonical SMILES after pH correction.")
    if Chem.MolFromSmiles(protomer) is None:
        raise ValueError(f"pH-corrected protomer is not parseable by RDKit: {protomer}")
    return protomer, f"OpenBabel pH-aware protonation (pH {ph:g}) + RDKit 3D embedding"


def prepare_ligand_from_smiles(
    smiles: str,
    ligand_name: str,
    output_dir: Path,
    ph: float = 7.4,
    require_protonation: bool = True,
    pubchem_cid: int | None = None,
) -> tuple[Path, Path, dict]:
    """Prepare a ligand SDF/PDBQT pair with explicit pH protonation.

    Returns ``(sdf_path, pdbqt_path, report)``.  The report records the input
    and prepared SMILES, formal charge, protonation method and the explicit
    limitations (no tautomer enumeration).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    smiles = (smiles or "").strip()
    if not smiles:
        raise ValueError("Empty SMILES string provided.")

    input_mol = Chem.MolFromSmiles(smiles)
    if input_mol is None:
        raise ValueError(f"Invalid SMILES string: '{smiles}'")

    prepared_smiles = smiles
    protonation_applied = False
    protonation_error = ""
    protonation_method = "none (RDKit AddHs only; pH not applied)"

    if require_protonation:
        try:
            prepared_smiles, protonation_method = _apply_ph_protonation(smiles, ph)
            protonation_applied = True
        except Exception as exc:
            protonation_error = str(exc)

    if require_protonation and not protonation_applied:
        raise ValueError(
            "Ligand protonation at pH 7.4 could not be applied, so the ligand was not prepared. "
            f"Reason: {protonation_error}. Install a working OpenBabel build or pass "
            "require_protonation=False only for an explicitly non-physiological experiment."
        )

    mol = Chem.MolFromSmiles(prepared_smiles)
    if mol is None:
        raise ValueError(f"Could not build the prepared ligand from SMILES: '{prepared_smiles}'")
    mol = _build_3d_conformer(mol)

    sdf_path = output_dir / f"{ligand_name}_3d.sdf"
    writer = Chem.SDWriter(str(sdf_path.resolve()))
    writer.write(mol)
    writer.close()

    pdbqt_path = output_dir / f"{ligand_name}.pdbqt"
    ob_mol = read_file(sdf_path, "sdf")
    write_file(ob_mol, pdbqt_path, "pdbqt")

    # Sanitize ligand PDBQT for any accidental invalid AutoDock atom types (e.g. Cm -> C)
    # and guarantee non-empty columns 78-79 with line length >= 79 to prevent parse_pdbqt.cpp(69)
    lines = []
    with open(pdbqt_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                clean_line = line.rstrip("\r\n")
                if len(clean_line) < 77:
                    clean_line = clean_line.ljust(77)
                ad_type = clean_line[77:79].strip() if len(clean_line) >= 79 else clean_line[77:].strip()
                if ad_type == "Cm":
                    ad_type = "C"
                if not ad_type:
                    atom_name = clean_line[12:16].strip() if len(clean_line) >= 16 else ""
                    alpha = "".join(ch for ch in atom_name if ch.isalpha())
                    ad_type = alpha[0].upper() if alpha else "C"
                type_formatted = ad_type.ljust(2) if len(ad_type) == 1 else ad_type[:2]
                clean_line = clean_line[:77].ljust(77) + type_formatted + "\n"
                lines.append(clean_line)
            else:
                lines.append(line)
    with open(pdbqt_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    num_rotatable_bonds = AllChem.CalcNumRotatableBonds(mol)
    formal_charge = int(Chem.GetFormalCharge(mol))
    canonical_prepared = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)

    report = {
        "ligand_name": ligand_name,
        "smiles": smiles,
        "prepared_smiles": canonical_prepared,
        "sdf_file": str(sdf_path.name),
        "pdbqt_file": str(pdbqt_path.name),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "total_atom_count": mol.GetNumAtoms(),
        "rotatable_bonds": num_rotatable_bonds,
        "formal_charge": formal_charge,
        "ph": ph,
        "pubchem_cid": pubchem_cid,
        "pubchem_cid_recorded": bool(pubchem_cid),
        "protonation_applied": protonation_applied,
        "protonation_method": protonation_method,
        "tautomers_enumerated": False,
        "protonation_note": (
            f"Protomer prepared at pH {ph:g} with {protonation_method}. "
            "Tautomers were not enumerated; the recorded prepared SMILES is the exact docked species."
            if protonation_applied
            else "pH protonation was not applied; this is not a physiological-pH preparation."
        ),
    }
    return sdf_path, pdbqt_path, report
