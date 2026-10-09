# Consolidated Repairs: Cross-Protein, Reference-Ligand Validation, and Session-State Isolation

Apply only the repairs in this document. Do not hard-code hemoglobin or any other protein/ligand as a required default. The application must remain general for any user-supplied structures.

## 1. Preserve provenance for Exploratory Cross-Protein Docking

**File:** `docking/comparison.py` — `run_exploratory_cross_docking()`

The cross-protein function creates run directories but must save complete provenance for each receptor exactly as the matched workflow does.

For each cross-protein receptor:

1. Save `receptor_source.pdb` immediately.
2. Write `manifest.initial.json` before receptor/ligand preparation or Vina execution.
3. On any preparation, Vina, or clustering failure, finalize both manifests with `status="run failed"`, exact failure reason, input fingerprint, receptor-specific grid, pocket list, and completed metadata.
4. On success, write `manifest.json` with the receptor's own preparation report, ligand report, protonation note, Vina settings, grid, clusters, per-seed output, and selected pose provenance.

Use the following status for successful cross-protein runs:

```python
status = "exploratory—cross-protein docking"
```

Do not store or display `delta_score`, a ranking, or a mutation-effect conclusion for this workflow.

Each final manifest must retain:

```python
{
    "input_fingerprint": input_fingerprint,
    "target": receptor_target_meta,
    "binding_site": {
        "method": "residues",
        "residues": receptor_pocket_residues,
        "padding": padding,
    },
    "ligand": {
        **ligand_preparation_report,
        "protonation_note": PROTONATION_ASSUMPTION_NOTE,
    },
    "results": {
        "top_score": min(p.score for p in poses),
        "dominant_cluster_id": selected_cluster.cluster_id,
        "dominant_cluster_median": selected_cluster.median_score,
        "dominant_cluster_seed_count": selected_cluster.seed_count,
        "clusters": serialized_clusters,
        "per_seed": per_seed_metadata,
    },
}
```

## 2. Do not call a one-seed pose a consensus pose

**Files:** `docking/comparison.py`, `docking/ui.py`

Add:

```python
MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT = 2
```

For both receptors in `run_exploratory_cross_docking()`:

```python
selected_cluster = _select_consensus_cluster(clusters)
cluster_reproducible = (
    selected_cluster.seed_count >= MIN_EXPLORATORY_CLUSTER_SEED_SUPPORT
)
```

Return and save:

```python
"cluster_reproducible": cluster_reproducible,
"dominant_cluster_seed_count": selected_cluster.seed_count,
```

UI behavior:

- If `cluster_reproducible` is true, write: `Reproducible exploratory cluster (>= 2 independent seeds)`.
- If false, write: `Exploratory pose without independent-seed support`.
- Do not call it a consensus pose and do not emphasize its median score as a stable result when support is below two seeds.

The 3D visualization and contact table must always use the same returned `dominant_cluster_pose` whose ID, seed, rank, score, and seed support are displayed. Never use `clusters[0]` as a shortcut.

## 3. Validate that the WT reference ligand is the same molecule as the submitted ligand

**Files:** `docking/comparison.py`, `docking/validation.py`, `docking/ui.py`

Before assigning:

```text
Protocol validated on WT reference
```

require a strict ligand-identity gate.

1. Obtain a verified reference SMILES for the selected co-crystal ligand from RCSB Chemical Component data, or require a verified reference SMILES field from the user.
2. Parse both the verified reference SMILES and the submitted SMILES using RDKit.
3. Compare all of the following after a documented standardization step:

```python
reference_canonical_isomeric_smiles == submitted_canonical_isomeric_smiles
reference_heavy_atom_count == submitted_heavy_atom_count
reference_molecular_formula == submitted_molecular_formula
```

4. If any check fails, do not calculate validation RMSD. Set:

```text
Not validated for affinity comparison: the submitted ligand is not the verified WT co-crystal ligand.
```

5. Save the reference component ID, reference SMILES SHA-256, submitted SMILES SHA-256, identity-gate result, and reason in the two manifests.

Only after the identity gate passes, run the existing no-fit, symmetry-aware RMSD validation. Make its MCS chemically strict:

```python
mcs = rdFMCS.FindMCS(
    [mol_ref, mol_pose],
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareOrderExact,
    ringMatchesRingOnly=True,
    completeRingsOnly=True,
    matchValences=True,
)
```

Keep the validation details in the manifest: RMSD method, threshold, Top-1 RMSD, Top-N recovery RMSD, supporting seeds, and pass/fail state.

## 4. Make selected chains authoritative in Cross-Protein mode

**File:** `docking/ui.py`

Cross-Protein mode lets the user select a chain for each receptor. Ensure entered pocket residues belong to the selected chain.

Add:

```python
def validate_pocket_chain(residue_specs: list[str], selected_chain: str, label: str) -> None:
    for spec in residue_specs:
        if ":" not in spec:
            raise ValueError(
                f"{label} pocket residue {spec!r} must use CHAIN:RESIDUE notation."
            )
        chain, _ = spec.split(":", 1)
        if chain.strip() != selected_chain:
            raise ValueError(
                f"{label} pocket residue {spec!r} does not belong to selected chain {selected_chain!r}."
            )
```

Call it for both receptor pocket lists before starting Vina. Save selected-chain IDs in manifest target metadata.

## 5. Clear stale protein data and all dependent results from `session_state`

**File:** `analysis/ui.py`

Loading a different PDB correctly changes docking fingerprints later, but clearing an input field or removing an uploaded file does not currently remove the old loaded protein state. Add a central invalidation helper:

```python
def invalidate_dependent_results(reason: str) -> None:
    """Clear displayed/in-memory protein-derived results, never saved run folders."""
    keys = {
        "h_results", "m_results", "h_selected_chain", "m_selected_chain",
        "last_coordinate_mode",
        "current_docking_result", "dock_single_input_fingerprint",
        "matched_comparison_result", "dock_comp_input_fingerprint",
        "exploratory_cross_result",
        "screening_results", "dock_screen_input_fingerprint",
        "comp_h_chain_sel", "comp_m_chain_sel",
        "cross_r1_chain_sel", "cross_r2_chain_sel",
        "comp_pocket_input", "comp_pocket_val",
        "comp_override_m_pocket_input", "comp_override_reason_input",
        "cross_r1_pocket_input", "cross_r2_pocket_input",
        "chk_comp_override_pocket",
        "comp_ref_lig_resname", "comp_ref_lig_chain", "comp_ref_lig_seq",
    }
    for key in keys:
        st.session_state.pop(key, None)
    st.session_state["last_result_invalidation_reason"] = reason
```

Call this helper after successfully replacing coordinates from a PDB ID or changed uploaded file, before `st.rerun()`.

For uploads, store a SHA-256 hash of `file.getvalue()` and invalidate only when that hash changes. This prevents an unchanged uploaded file from repeatedly clearing results during normal reruns.

Add explicit `Clear Healthy protein` and `Clear Mutant protein` buttons:

```python
def clear_protein(prefix: str, label: str) -> None:
    for key in {
        f"{prefix}_pdb", f"{prefix}_assembly_pdb", f"{prefix}_id",
        f"{prefix}_id_in_val", f"{prefix}_id_widget",
        f"{prefix}_upload_sha256", f"{prefix}_selected_chain",
        f"{prefix}_results",
    }:
        st.session_state.pop(key, None)
    invalidate_dependent_results(f"{label} protein was cleared by the user.")
```

After clearing, call `st.rerun()`. The Docking tab must show its normal “load both proteins first” message and must not display an old run.

Never delete `docking_runs/` or manifest files during invalidation. Those old folders are separate provenance records and must remain available.

Show one confirmation message after invalidation:

```python
reason = st.session_state.pop("last_result_invalidation_reason", None)
if reason:
    st.info(
        f"Previous analysis and docking results were cleared to prevent mixing datasets. "
        f"{reason} Saved run folders remain available as separate provenance records."
    )
```

Keep the existing coordinate-hash fingerprints in `docking/ui.py` as a second safety layer.

## Required regression tests

Add tests for all of the following:

1. Successful mocked Cross-Protein docking creates two final manifests with their own pocket lists and no `delta_score`.
2. Failure during the second Cross-Protein Vina run finalizes both manifests as `run failed` with the exact failure reason.
3. A one-seed Cross-Protein cluster is labelled non-reproducible; a two-seed cluster is labelled reproducible.
4. Two ligands with equal heavy-atom counts but different connectivity or bond order fail the strict reference-ligand identity gate and cannot receive `Protocol validated on WT reference`.
5. Cross-Protein pocket residues belonging to a chain other than the selected chain are rejected.
6. Replacing only healthy coordinates clears `current_docking_result`, `matched_comparison_result`, `exploratory_cross_result`, and `screening_results`.
7. Replacing only mutant coordinates clears the same result objects.
8. Clearing either protein removes its source PDB/chain and makes both-protein docking unavailable.
9. Re-uploading identical bytes does not invalidate current results; changed bytes do.
10. Old manifest files remain on disk after session-state invalidation.

## Completion commands

Run from the `Bio_V2` project root:

```powershell
& "C:\Users\hassan\AppData\Local\Programs\Python\Python312\python.exe" -m compileall -q b.py analysis docking tests
& "C:\Users\hassan\AppData\Local\Programs\Python\Python312\python.exe" -m unittest discover -s tests -v
& "C:\Users\hassan\AppData\Local\Programs\Python\Python312\python.exe" -c "from streamlit.testing.v1 import AppTest; at=AppTest.from_file('b.py'); at.run(timeout=45); assert not at.exception, [str(x.value) for x in at.exception]; print('AppTest OK')"
```
