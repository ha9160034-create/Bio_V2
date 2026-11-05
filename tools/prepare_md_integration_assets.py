"""Prepare real 3HTB/JZ4 installation-test inputs using the existing Vina workflow.

Workflow sources: wiki/gromacs-protein-ligand-protocol.md,
wiki/autodock-vina.md, wiki/docking-failure-modes.md,
wiki/deep-learning-docking.md. A single-seed smoke test does not validate
pose recovery, binding affinity, or MD convergence.
"""

import hashlib
import json
import sys
from dataclasses import asdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from docking.binding_site import define_grid_from_ligand
from docking.engine import run_vina_single_seed
from docking.ligand_prep import prepare_ligand_from_smiles
from docking.rcsb import fetch_chem_comp_smiles, fetch_deposited_pdb_text
from docking.receptor_prep import prepare_receptor


def main():
    base = PROJECT_ROOT / "tmp_gmx_test"
    run = base / "integration"
    marker = run / "asset_setup.json"
    if run.exists() and any(run.iterdir()) and not marker.is_file():
        raise RuntimeError(f"Refusing to overwrite pre-existing integration assets: {run}")
    if marker.is_file():
        previous = json.loads(marker.read_text(encoding="utf-8"))
        if previous.get("purpose") != "real_3htb_jz4_installation_test":
            raise RuntimeError("Integration asset ownership marker does not match.")
        if previous.get("status") == "completed":
            required = ["manifest.json", "receptor_clean.pdb", "out_seed_42.pdbqt", "JZ4_3d.sdf"]
            if all((run / name).is_file() for name in required):
                print(f"Real integration assets already prepared: {run}")
                return

    run.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps({"purpose": "real_3htb_jz4_installation_test", "status": "preparing"}), encoding="utf-8")
    pdb_path = base / "3HTB.pdb"
    if pdb_path.exists():
        pdb_text = pdb_path.read_text(encoding="utf-8")
        pdb_source = "existing local reference file"
    else:
        pdb_text = fetch_deposited_pdb_text("3HTB")
        if not pdb_text:
            raise RuntimeError("Official 3HTB coordinates could not be downloaded.")
        pdb_path.write_text(pdb_text, encoding="utf-8")
        pdb_source = "https://files.rcsb.org/download/3HTB.pdb"
    smiles = fetch_chem_comp_smiles("JZ4")
    if not smiles:
        raise RuntimeError("The official JZ4 chemical-component SMILES is unavailable.")
    instances = sorted({
        (line[21:22], line[22:26].strip())
        for line in pdb_text.splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == "JZ4"
    })
    if not instances:
        raise RuntimeError("The downloaded structure has no JZ4 reference ligand.")
    chain, resseq = instances[0]
    grid, _coords, _reference = define_grid_from_ligand(pdb_text, "JZ4", chain_id=chain, resseq=resseq)
    receptor, receptor_pdbqt, receptor_report = prepare_receptor(pdb_text, run, reference_ligand_resname="JZ4")
    sdf, ligand_pdbqt, ligand_report = prepare_ligand_from_smiles(smiles, "JZ4", run)
    poses, vina_report = run_vina_single_seed(receptor_pdbqt, ligand_pdbqt, grid, run)
    if not poses:
        raise RuntimeError("Vina produced no actual poses for the installation test.")
    manifest = {
        "run_id": "integration_3htb_jz4",
        "status": "completed",
        "purpose": "installation_smoke_test_not_affinity_prediction",
        "target": {"pdb_id": "3HTB", "coordinate_mode": "deposited"},
        "ligand": {**ligand_report, "name": "JZ4"},
        "grid": asdict(grid),
        "grid_reference_instance": {"chain": chain, "resseq": resseq},
        "receptor_preparation": receptor_report,
        "vina_runs": [vina_report],
        "protocol_validation": {"status": "not_validated", "reason": "single_seed_installation_smoke_test"},
        "sources": {"pdb": pdb_source, "ligand": "https://data.rcsb.org/rest/v1/core/chemcomp/JZ4"},
        "sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in (pdb_path, receptor, sdf)},
    }
    (run / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    marker.write_text(json.dumps({"purpose": "real_3htb_jz4_installation_test", "status": "completed"}), encoding="utf-8")
    print(f"Real Vina smoke test completed: {len(poses)} poses; inputs saved in {run}")


if __name__ == "__main__":
    main()
