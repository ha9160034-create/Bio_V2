#!/usr/bin/env python3
"""Headless matched comparison: 1M17 (WT EGFR) vs 2JIT (T790M), erlotinib.

Produces the same provenance manifests as the Streamlit comparison tab and a
compact JSON summary for review.  Run from the repository root:

    python scripts/run_t790m_comparison.py

The mutant arm is apo (2JIT), so its manifest is expected to carry
"matched comparison-mutant arm (WT validation not transferable)".
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path

from docking.comparison import run_matched_docking_comparison
from docking.pair_validation import validate_matched_pair
from docking.presets import get_preset
from docking.rcsb import fetch_deposited_pdb_text, fetch_entry_quality

ROOT = Path(__file__).resolve().parents[1]
HEALTHY_RUN = ROOT / "docking_runs" / "dock_1M17_20261115_164801_5913a8"


def _pocket_list(spec: str) -> list[str]:
    return [item.strip() for item in spec.split(",") if item.strip()]


def main() -> int:
    preset = get_preset("egfr_erlotinib_t790m")
    healthy_text = (HEALTHY_RUN / "receptor_source.pdb").read_text(
        encoding="utf-8", errors="replace"
    )
    mutant_text = fetch_deposited_pdb_text(preset["mutant_id"])
    if not mutant_text:
        raise SystemExit("Could not download 2JIT from RCSB.")

    pair = validate_matched_pair(
        healthy_text,
        mutant_text,
        preset["healthy_chain"],
        preset["mutant_chain"],
    )
    print(
        f"Pair identity: {pair['identity_pct']:.2f}% | "
        f"valid={pair['is_valid']} | mutations={pair['mutation_count']}"
    )
    if not pair["is_valid"]:
        raise SystemExit(f"Pair validation failed: {pair.get('error')}")

    quality = fetch_entry_quality(preset["mutant_id"]) or {}
    healthy_meta = {
        "pdb_id": preset["healthy_id"],
        "chain": preset["healthy_chain"],
        "title": "EGFR kinase domain with erlotinib (AQ4)",
        "source": "RCSB PDB",
        "resolution_angstrom": 2.6,
        "r_free": 0.295,
    }
    mutant_meta = {
        "pdb_id": preset["mutant_id"],
        "chain": preset["mutant_chain"],
        "title": "EGFR T790M kinase domain (apo)",
        "source": "RCSB PDB",
        "resolution_angstrom": quality.get("resolution_angstrom"),
        "r_free": quality.get("r_free"),
    }

    pocket = _pocket_list(preset["pocket_residues"])
    mutant_pocket = []
    for spec in pocket:
        mapped = pair["mapping"].get(spec) or pair["mapping"].get(spec.split(":")[-1])
        if not mapped:
            raise SystemExit(f"Could not map healthy pocket residue {spec} to 2JIT.")
        mutant_pocket.append(mapped)
    print(f"Healthy pocket ({preset['healthy_id']}): {pocket}")
    print(f"Mapped mutant pocket ({preset['mutant_id']}): {mutant_pocket}")
    ref_pdb = (HEALTHY_RUN / "reference_ligand.pdb").read_text(
        encoding="utf-8", errors="replace"
    )

    comparison = run_matched_docking_comparison(
        healthy_pdb_text=healthy_text,
        mutant_pdb_text=mutant_text,
        healthy_target_meta=healthy_meta,
        mutant_target_meta=mutant_meta,
        ligand_smiles=preset["drug_smiles"],
        ligand_name="Erlotinib",
        healthy_pocket_residues=pocket,
        mutant_pocket_residues=mutant_pocket,
        target_domain=preset["target_domain"],
        padding=8.0,
        seeds=[42, 101, 2024],
        exhaustiveness=8,
        pair_validation_meta=pair,
        wt_reference_ligand_meta={
            "ref_ligand_pdb_text": ref_pdb,
            "resname": "AQ4",
            "ref_smiles": preset["drug_smiles"],
        },
    )

    summary = {
        "generated_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "pair_validation": {
            key: value
            for key, value in pair.items()
            if key not in {"mapping", "mutations"}
        },
        "delta_score": comparison["delta_score"],
        "delta_score_uncertainty_kcal": comparison[
            "delta_score_uncertainty_kcal"
        ],
        "validation_status": comparison["validation_status"],
        "validation_reason": comparison["validation_reason"],
        "binding_mode_status": comparison["binding_mode_status"],
        "mutant_validation_status": comparison["mutant_validation_status"],
        "mutant_validation_reason": comparison["mutant_validation_reason"],
        "binding_mode_check": comparison["binding_mode_check"],
        "healthy": {
            key: comparison["healthy"][key]
            for key in (
                "run_id",
                "target_id",
                "dominant_cluster_id",
                "dominant_cluster_top_score",
                "dominant_cluster_median",
                "dominant_cluster_seeds",
                "seed_count",
            )
        },
        "mutant": {
            key: comparison["mutant"][key]
            for key in (
                "run_id",
                "target_id",
                "dominant_cluster_id",
                "dominant_cluster_top_score",
                "dominant_cluster_median",
                "dominant_cluster_seeds",
                "seed_count",
            )
        },
        "warnings": comparison["warnings"],
        "disclaimer": comparison["disclaimer"],
    }
    out = ROOT / "docking_runs" / "comparison_t790m_20261118.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
