# Generic Matched Docking: Scientific-Reliability Repairs

## Scope and non-negotiable design rule

The application must remain usable with **any user-supplied pair of structures**. Do not hard-code hemoglobin, 2HHB, 2HBS, Voxelotor, Erlotinib, `A:790`, `A:858`, or any other protein-specific values as requirements or universal defaults.

Hemoglobin is only a regression example. It must never be a condition for running the UI or for accepting another protein pair.

The current 2HHB/2HBS calculation completed correctly as a Vina calculation: the dominant cluster median was `-5.79 kcal/mol` for each structure and the displayed `ΔScore = 0.00 kcal/mol` follows the program's formula. It is still correctly labelled `exploratory—no reference validation`; it is **not** evidence that a drug has no biological effect.

## 1. Separate “matched variant comparison” from “exploratory cross-protein docking”

**Files:** `docking/ui.py`, `docking/comparison.py`, optionally a new `docking/pair_validation.py`

At present, a free-text `target_domain` and identical ligand/settings are the only biological checks. The code can therefore call two unrelated proteins a matched healthy/mutant comparison if both happen to contain the entered residue labels. That makes `ΔScore` scientifically misleading.

### Required behavior

Offer two explicit workflows:

1. **Matched WT-versus-variant comparison**
   - Only for the same protein/gene/domain, with the two selected chains homologous.
   - Produce a protocol-specific `ΔScore` only after pair validation passes.
2. **Exploratory cross-protein docking**
   - Allows any two structures and independently specified pockets.
   - Display two separate docking summaries.
   - Never display a matched `ΔScore`, “more/less favorable for the mutant”, ranking, or drug-response conclusion.

### Matched-pair validation

Use the chains selected in the Structural Analysis tab (`h_selected_chain` and `m_selected_chain`) and their active coordinate texts. Before invoking Vina:

1. Extract their protein sequences.
2. Perform global sequence alignment with the existing `analysis.alignment.get_alignment` implementation.
3. Calculate identity only over aligned, non-gap positions.
4. Require at least `90.0%` identity by default. Put the threshold in a named constant such as `MIN_MATCHED_PAIR_IDENTITY = 90.0` and display it.
5. Build an explicit healthy-to-mutant residue mapping from the alignment, including insertion codes.
6. Verify every selected healthy pocket residue maps to a real mutant residue in the selected mutant chain. Do not assume the same chain ID or residue number in both structures.
7. Save the selected chains, identity, aligned-residue count, mutation count, and residue mapping in both manifests.

If any check fails, do not run the matched calculation. Show a clear message such as:

```text
Matched comparison was not run: selected chains are not a validated WT/variant pair
(sequence identity 34.6%, required >= 90.0%). Use Exploratory Cross-Protein Docking instead.
```

This preserves the user's freedom to dock any proteins while preventing invalid biological comparisons from being reported as mutant effects.

## 2. Remove unsafe protein-specific defaults from Matched Comparison

**File:** `docking/ui.py`, `_render_matched_comparison_ui`

The current matched workflow falls back to `A:790, A:858`, `Kinase Domain`, and Erlotinib. These belong to one EGFR example and are unsafe for a general pair such as hemoglobin.

### Required change

Use blank values unless the user has explicitly entered a value in this session:

```python
default_pocket = st.session_state.get("comp_pocket_val", "")
default_domain = st.session_state.get("comp_domain_val", "")
default_l_name = st.session_state.get("comp_lig_name_val", "")
default_smiles = st.session_state.get("comp_smiles_val", "")
```

Require all four fields before a matched run. Do not copy `dock_pocket_res_val`, `dock_lig_name_val`, or `dock_lig_smiles_val` automatically from the single-docking form, because that can silently reuse a ligand or pocket intended for another target.

An optional **Load documented example** button may populate a named preset after user action. It must state the PDB IDs, ligand, pocket source, and that it is an example only.

## 3. Let each receptor use its own mapped pocket definition

**Files:** `docking/ui.py`, `docking/comparison.py`

The current implementation sends one `pocket_residues` list to both receptors. This fails or becomes scientifically wrong when equivalent chains use different chain IDs, residue numbering, or insertion codes.

### Required API change

Change the comparison function to accept two validated lists:

```python
healthy_pocket_residues: list[str]
mutant_pocket_residues: list[str]
```

Build each grid from its own list:

```python
h_grid = define_grid_from_residues(
    healthy_pdb_text, healthy_pocket_residues, padding=padding
)
m_grid = define_grid_from_residues(
    mutant_pdb_text, mutant_pocket_residues, padding=padding
)
```

Store both lists and the residue mapping in the manifest. In the UI, let the user enter the WT residues and show the mapped mutant residues as a reviewable table. Allow manual correction only with an explicit reason recorded in the manifest.

## 4. Use the same consensus pose for scores, visualization, and contacts

**File:** `docking/comparison.py` and `docking/ui.py`, around current lines 1322-1323.

`verify_comparison_clusters()` selects the reproducible consensus cluster by seed support, then median score. The UI currently calculates contact candidates from `clusters[0]`, which may be a different cluster. The displayed contacts/3D pose can therefore disagree with the `ΔScore` shown above them.

### Required change

Return the selected cluster identity and representative pose source explicitly:

```python
"dominant_cluster_id": h_top_cluster.cluster_id,
"dominant_cluster_pose": h_top_cluster.representative_pose,
```

and the equivalent mutant values. Then use only those returned values in the UI:

```python
h_top_pose = comp["healthy"]["dominant_cluster_pose"].pdbqt_block
m_top_pose = comp["mutant"]["dominant_cluster_pose"].pdbqt_block
```

Label the rendered pose with its cluster ID, seed, rank, seed support, and median score. Do not use `clusters[0]` for any result claimed to describe the consensus comparison.

## 5. Make validation status scientifically meaningful for a matched comparison

**Files:** `docking/ui.py`, `docking/comparison.py`, `docking/storage.py`

Three seeds establish only **computational reproducibility**. They do not establish that the predicted binding pose or score is accurate.

For a Matched WT-versus-variant result, show these separate states:

| Condition | Required displayed status |
|---|---|
| No suitable co-crystal reference for the WT protocol | `Exploratory: no pose validation` |
| WT co-crystal redocking succeeds (Top-1 or documented Top-N RMSD threshold, with multi-seed support) | `Protocol validated on WT reference` |
| WT redocking fails or is unavailable | `Not validated for affinity comparison` |
| Selected consensus cluster has fewer than two supporting seeds | Reject `ΔScore` |

Add optional WT reference-ligand selection (resname, chain, residue number) to Matched Comparison. When supplied, redock it using the exact receptor preparation, grid, ligand preparation, Vina version, exhaustiveness, seeds, number of modes, and energy range used for the comparison. Compute symmetry-aware heavy-atom RMSD without fitting, using the existing validation utilities.

Only a passed WT redocking validates this **protocol family**. It does not turn Vina score into experimental affinity or clinical evidence.

## 6. Correctly communicate what every displayed number means

### Values that are computationally correct in the current code

- **Vina score:** correctly read from Vina output for the stated prepared receptor, ligand, grid, seed, and search settings.
- **Dominant-cluster median:** a reproducibility-oriented summary, because the selected cluster prioritizes independent seed support before its median Vina score.
- **ΔScore:** mutant dominant-cluster median minus healthy dominant-cluster median.
- **RMSD:** meaningful only when the docked and reference ligands are chemically identical and a mapped, no-fit heavy-atom RMSD calculation is available.

### Values that must never be claimed

- Experimental binding free energy (`ΔG`)
- Drug efficacy, resistance, dose, safety, or clinical benefit
- A genuine hydrogen bond merely from the distance-based contact candidates
- A cross-protein or cross-pocket score difference as though it were a WT-versus-mutant effect

The current ligand-preparation parameter `ph=7.4` is metadata only; it does not enumerate or prove the dominant protonation/tautomer state. Add this exact note to the UI and manifest:

```text
The submitted SMILES was docked as prepared. “pH 7.4” is an assumption label,
not an ionization or tautomer prediction.
```

For metalloproteins/cofactor systems, retain the existing warning: standard Gasteiger charges are not specialized metal/cofactor parameterization. Results must be labelled exploratory unless independently validated.

## 7. Required regression tests

Add tests that are generic; none may require hemoglobin as a UI default.

1. Two unrelated selected chains must be rejected in Matched mode and offered the exploratory cross-protein workflow.
2. A high-identity synthetic WT/variant pair with a residue-number shift must create correct separate WT and mutant pocket lists.
3. The UI contact pose must be the returned consensus representative, not `clusters[0]` when they differ.
4. A cluster supported by one seed must prevent `ΔScore` output.
5. A matched run without reference redocking must be labelled `Exploratory: no pose validation`.
6. A passing synthetic or documented redocking case must show its RMSD criterion, Top-1/Top-N result, seed support, and no-fit mapping method.
7. Empty matched fields must block execution; no implicit EGFR, Erlotinib, or residue defaults may be used.
8. The existing 2HHB/2HBS case may remain an optional regression fixture for the historical `Cm` parser fault, but it must never be a required runtime protein, default UI input, or universal validation benchmark.

## Required completion checks

```powershell
& "C:\Users\hassan\AppData\Local\Programs\Python\Python312\python.exe" -m compileall -q b.py analysis docking tests
& "C:\Users\hassan\AppData\Local\Programs\Python\Python312\python.exe" -m unittest discover -s tests -v
& "C:\Users\hassan\AppData\Local\Programs\Python\Python312\python.exe" -c "from streamlit.testing.v1 import AppTest; at=AppTest.from_file('b.py'); at.run(timeout=30); assert not at.exception, [str(x.value) for x in at.exception]; print('AppTest OK')"
```

Do not report the workflow as scientifically validated merely because these software tests pass. Software tests prove implementation behavior; redocking and appropriate external evidence validate a docking protocol.
