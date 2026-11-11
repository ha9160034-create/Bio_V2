# Changelog

All notable changes to the BIOBIO runtime and packaging. Format follows
Keep a Changelog; versions follow SemVer.

## [1.0.0-rc2] — 2026-11-10

Phase 2 of `PLAN.md`: deterministic ligand topology and charge neutrality.

### Added
- `_ligand_charge_sum()` + `_enforce_ligand_charge()` in `md/topology.py`:
  the production ITP is rejected with `TopologyBuildError` when
  `|sum(q) − net_charge| > 1e-4 e`, before position restraints or `grompp`.
  `parameterize_ligand_acpype()` now returns `charge_sum`, `charge_atoms` and
  `net_charge` for traceability.
- Explicit exit codes in `acpype_wrapper.py`: `2` usage/environment,
  `3` ACPYPE failed without Amber files, `4` no usable output, `5` charge
  validation, `6` ITP/GRO inconsistency, `7` converter failure.

### Changed
- `acpype_wrapper.py` rewritten around the file-based decision from
  `PLAN.md` §2.1: ACPYPE's exit code is logged but never discards physically
  valid output; a valid GMX pair is never overwritten; complete Amber files
  are converted with ParmEd when the GMX pair is missing, inconsistent or
  non-neutral.
- `convert_amber_to_gmx.py`: charge balancing now mirrors
  `acpype.topol.Topology.balanceCharges()` exactly — the residual is spread
  over **every** atom holding the maximum `|q|`, not a single host atom.
  `--net-charge` defaults to `round(sum(q))` (ACPYPE's own rule) instead of 0.
- ITP/GRO integrity is checked before acceptance: non-empty files, a valid
  `[ atoms ]` block, and equal ITP/GRO atom counts.
- Removed the `acpype/utils.py` monkey-patch (`PLAN.md` P1-7). `md/topology.py`
  now exports `AMBERHOME`/`PATH` from the runtime and leaves ACPYPE unmodified.
- `environment.lock` header refreshed for the rc2 image; the re-exported conda
  solution (214 URLs) and `pip freeze` are identical to rc1.

### Verification
- Real ACPYPE run on `ERLOT_pose.mdl` (AM1-BCC, GAFF2, 1m40s):
  `atoms=52`, `sum(q)=+0.000000 e`; wrapper accepted the output unchanged.
- Re-ran ACPYPE after removing the P1-7 patch: charges, atom types and
  coordinates are identical; the only file difference is the timestamp comment.
- Baseline Amber files (`ERLOT_AC.prmtop`): `residual_before=+0.00200000 e`
  → `sum_after=+0.00000000 e` via the ParmEd fallback.
- Negative test: a +1 e ligand with `net_charge=0` raises
  `TopologyBuildError` (`NEGATIVE_GUARD_PASSED`).
- Wrapper state tests (fake ACPYPE, real converter):
  `rc=1 + Amber → 0`, `rc=1 no Amber → 3`, `rc=0 no output → 4`,
  `rc=0 valid → 0 (ITP hash unchanged)`, `+1 e GMX no Amber → 5`.
- `grompp` on the archived 102,128-atom system with the new ITP:
  **0 warnings** (baseline had 1: PME net charge), tpr produced. Repeated with
  the no-patch ACPYPE output: same result.
- `verify_stack.py`: 25/25; `import md.simulation, md.topology, docking.storage`
  succeeds.
- Image: `biobio-app:1.0.0-rc2` = `sha256:43002bd7…`, 5.80 GB (local only).

## [1.0.0-rc1] — 2026-11-10

Phase 1 of `PLAN.md`: reproducibility of the scientific runtime.

### Added
- `constraints.txt` — single source of truth for every pip-installed version.
- `environment.lock` — exact conda + pip solution exported from the built image.
- `expected_versions.json` + `scripts/verify_stack.py` — fail-closed 25-check gate.
- `scripts/export_environment.sh` — regenerates `environment.lock`.
- Vina SHA-256 enforcement in `docking/storage.py`
  (`fa0126a2…` on Linux; a mismatch raises `VinaBinaryError` and refuses to dock).
- Vina provenance in every docking manifest: `vina_path`, `vina_sha256`,
  `vina_sha256_expected`, `vina_sha256_match`.
- Image digest sidecar support in `Launcher.cs` (`biobio-image.digest`).

### Changed
- Base image pinned: `mambaorg/micromamba:2.10.0@sha256:e1476182…c1198`
  (replaces `:latest`; measured base is 2.10.0 / Debian 13).
- Conda pinned: `python=3.11.18`, `gromacs=2026.3=nompi_h26635d9_100`,
  `parmed=4.3.1=py311h74232ed_1`, `ambertools=26.0=cuda_None_nompi_py311h13e84fe_101`,
  `acpype=2023.10.27=pyhd8ed1ab_0`.
- Removed `autodock-vina` from apt: it provided a second Vina (1.2.7) that
  resolved first in `PATH`, while the app used 1.2.5.
- Removed conda `openbabel`; `openbabel-wheel==3.1.1.23` is now the only copy
  (conda and wheel share files, so keeping both made uninstall corrupt the other).
- `requirements.txt` is now top-level names only; versions live in `constraints.txt`.
- `Launcher.cs`: image tag `hassancreator/biobio-app:1.0.0` (no `latest`).

### Verification
- `verify_stack.py`: **25/25 checks passed**.
- Redocking 1M17 / Erlotinib, seed 42: scores `-6.822 / -6.656 / -6.622`,
  output PDBQT SHA-256 `1585e6e7…fa013` — **bit-identical to baseline**.
- EM from archived `em.tpr`: completes; 10-step run exits cleanly.
- Negative test: tampered Vina binary rejected with `VinaBinaryError`.
- Container boots, Streamlit answers HTTP 200, `import md.simulation` succeeds.
- Image size: 5.80 GB (baseline 6.09 GB).
