"""Stable ligand atom provenance.

References: wiki/ligand-preparation-and-protonation.md and
wiki/open-force-field-and-ligand-parameterization.md. Match coordinates only
between the undocked SDF/PDBQT pair; after Vina, use preserved serial numbers.
"""
import json

import numpy as np
from scipy.optimize import linear_sum_assignment

MAP_PROPERTY = "BIO_PDBQT_SERIAL_TO_TEMPLATE"


def record_atom_mapping(molecule, pdbqt_text: str) -> dict[str, int]:
    from rdkit import Chem

    heavy_mol = Chem.RemoveHs(molecule)
    template_xyz = heavy_mol.GetConformer().GetPositions()
    serials, pdbqt_xyz = [], []
    for line in pdbqt_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_type = line[77:].strip().upper()
        if atom_type.startswith("H"):
            continue
        serials.append(int(line[6:11]))
        pdbqt_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if len(pdbqt_xyz) != len(template_xyz) or len(set(serials)) != len(serials):
        raise ValueError("SDF/PDBQT atom mapping requires equal, unique heavy atoms.")
    distances = np.linalg.norm(np.asarray(pdbqt_xyz)[:, None, :] - template_xyz[None, :, :], axis=2)
    rows, columns = linear_sum_assignment(distances)
    if np.max(distances[rows, columns], initial=0) > 0.002:
        raise ValueError("Undocked SDF/PDBQT coordinates do not share the same frame.")
    mapping = {str(serials[row]): int(column) for row, column in zip(rows, columns)}
    molecule.SetProp(MAP_PROPERTY, json.dumps(mapping, sort_keys=True))
    return mapping
