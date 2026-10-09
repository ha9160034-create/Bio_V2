# Repair: Preserve Scientific Validation Status in Docking Manifests

## Problem

The matched-docking workflow can calculate a valid scientific outcome such as
`Protocol validated on WT reference`, but `docking/storage.py` does not accept
that text as a valid manifest status. The storage layer then silently converts
it to an exploratory status. This does not change a Vina score or RMSD, but it
corrupts the saved provenance and can make a validated protocol look
unvalidated after the page is refreshed or the run is opened later.

Do **not** edit or overwrite old `docking_runs` folders. They are historical
provenance records. Apply this repair only to new runs created after the code
change.

## Required implementation

### 1. Use fixed status values, not free-form explanatory sentences

In `docking/storage.py`, define one shared set of allowed manifest statuses.
Use exactly these normalized strings:

```python
MANIFEST_STATUSES = {
    "validation passed",
    "validation passed (Top-10 recovery)",
    "validation failed",
    "validation unavailable",
    "protocol validated on WT reference",
    "not validated for affinity comparison",
    "exploratory—no reference validation",
    "exploratory—insufficient independent seeds",
    "exploratory—cross-protein docking",
    "run failed",
}
```

Every call to `build_and_save_manifest(...)` must pass one of these values.
The status is a machine-readable scientific classification, not a long user
message.

### 2. Preserve the reason separately

In `docking/comparison.py`, replace variable-length values such as:

```python
"Not validated for affinity comparison: the submitted ligand is not ..."
```

with the fixed status:

```python
validation_status = "not validated for affinity comparison"
validation_reason = (
    "The submitted ligand is not chemically identical to the verified "
    "WT co-crystal reference ligand."
)
```

For a successful reference validation, use:

```python
validation_status = "protocol validated on WT reference"
validation_reason = (
    "The submitted ligand matched the verified WT reference ligand and "
    "the configured redocking validation criterion passed."
)
```

For no usable reference ligand, use:

```python
validation_status = "exploratory—no reference validation"
validation_reason = "No verified WT co-crystal ligand was available for redocking validation."
```

For a failed identity gate or failed RMSD/pose criterion, use:

```python
validation_status = "not validated for affinity comparison"
```

Do not use `run failed` for a completed calculation whose scientific validation
did not pass. `run failed` is reserved for an execution error.

Store the explanatory reason in the existing validation metadata, for example:

```python
results_meta["validation"] = {
    **results_meta.get("validation", {}),
    "status": validation_status,
    "reason": validation_reason,
}
```

Pass the same fixed `validation_status` into the final manifest. Return the
reason to the UI as a separate field and display it beneath the status.

### 3. Stop silent status coercion in `docking/storage.py`

Replace the current fallback that silently turns an unknown status into an
exploratory status (or `run failed`) with explicit normalization.

Keep only a deliberate legacy alias if needed by an old test or old caller:

```python
LEGACY_STATUS_ALIASES = {
    "exploratory": "exploratory—no reference validation",
}


def normalize_manifest_status(status: str) -> str:
    normalized = str(status).strip().lower()
    normalized = LEGACY_STATUS_ALIASES.get(normalized, normalized)
    if normalized not in MANIFEST_STATUSES:
        raise ValueError(f"Unsupported docking manifest status: {status!r}")
    return normalized
```

Then use `normalize_manifest_status(status)` before writing the manifest.
This makes a future typo visible during testing instead of silently changing
the scientific interpretation of a saved result.

If the code intentionally needs a new scientific state, add it to
`MANIFEST_STATUSES` first and add a matching test.

### 4. Keep cross-protein results unchanged

`run_exploratory_cross_docking(...)` must continue to save exactly:

```python
"exploratory—cross-protein docking"
```

It must not contain a mutation-effect `DeltaScore`, a WT/MT superiority claim,
or a `protocol validated on WT reference` label.

### 5. Correct one UI wording issue

The cross-protein workflow executes separate receptor calculations; it is not
currently parallel. Change any label reading:

```text
Independent parallel simulations
```

to:

```text
Independent separate simulations
```

Do not claim parallel execution unless the implementation actually dispatches
the calculations concurrently and safely collects both failures/results.

## Required tests

Add tests in `tests/test_repairs.py` (or the project’s docking test module):

1. Save a manifest with `protocol validated on WT reference`; assert that the
   JSON retains exactly that value.
2. Save a manifest with `not validated for affinity comparison`; assert that it
   remains that value and its `validation.reason` remains present.
3. Mock a successful matched WT-reference validation and assert that both the
   returned result and final manifest use `protocol validated on WT reference`.
4. Mock a reference-identity mismatch and assert that the final manifest is
   `not validated for affinity comparison`, never `run failed` and never an
   exploratory reference-validation status.
5. Assert that a cross-protein manifest remains
   `exploratory—cross-protein docking`.
6. Assert that an unknown status raises `ValueError` (or is handled by the one
   documented legacy alias only).

Run:

```powershell
& 'C:\Users\hassan\AppData\Local\Programs\Python\Python312\python.exe' -m compileall -q b.py analysis docking tests
& 'C:\Users\hassan\AppData\Local\Programs\Python\Python312\python.exe' -m unittest discover -s tests -v
```

The repair is accepted only if all tests pass and a newly generated matched
comparison manifest keeps its exact fixed scientific status after rerunning
the Streamlit page.
