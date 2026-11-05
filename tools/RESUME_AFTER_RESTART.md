# Bio_V2: resume after Windows restart

Updated: 2026-11-05 (Asia/Baghdad).

Project: `C:\Users\krar9\Desktop\Bio_V2-main`.

## Current state after the restart

- The user restarted Windows. WSL 2.7.14 and the running hypervisor were verified afterwards; do not ask for the same restart again.
- No pre-existing Linux distributions were found.
- The independent, limited-user Windows task `BioV2-MD-Backend-Setup` is running `tools/setup_md_backend.ps1` to install Ubuntu and then the MD runtime packages. It does not store credentials or restart Windows.
- The current authoritative installation status is `tools/md_backend_setup_status.json`, with details in `tools/md_backend_setup.log`. The older `wsl_setup_status.json` records the previous pre-reboot phase.
- A real AutoDock Vina run on public 3HTB/JZ4 reference inputs completed and produced 9 poses. Assets and their provenance are in `tmp_gmx_test/integration/`. This single-seed installation test does not establish binding affinity or pose-recovery validation.
- `tests.test_md_ui` passed using these assets and synthetic, explicitly labelled MD result fixtures. This is UI verification, not a real MD simulation.
- The Codex thread follow-up `bio-v2` is active again to continue verification after the background setup completes. Routine pending states should remain quiet.

## Saved and verified

- Code repairs are saved on disk: energy drift uses actual timestamps; missing analysis is not reported as zero; ligand RMSD selects heavy atoms; requested force-field versions are preserved; exported scripts use the custom index and configured salt.
- Tests in the project environment: 149 discovered, 142 passed, 7 skipped. These do not verify a real MD simulation.
- Project requirements match; pip check passed. The .venv reuses compatible system packages and overrides NumPy 2.2.6, Streamlit 1.62.0 and ACPYPE 2026.10.15. RDKit is pinned to 2026.3.6.
- Original files are backed up under `repair_backups/20261105-111508/`.
- Microsoft WSL 2.7.14 x64 MSI installed successfully (MSI exit 0).
- WSL and VirtualMachinePlatform Windows features enabled without an automatic restart.
- `tools/wsl_setup_status.json` reports `restart_required`.

## Next steps

1. Inspect the current setup task and status before starting any installer; do not duplicate a running download or setup process.
2. Preserve all existing distributions and data. Windows has already been restarted and WSL verified.
3. Allow the setup task to finish Ubuntu and GROMACS/ACPYPE/AmberTools installation in the same Linux backend; inspect and resolve any reported failure.
4. Verify backend discovery from the project's Python environment. Prepare actual reference assets if absent and run a short real MD integration test according to AGENTS.md and the relevant wiki protocol pages.
5. Launch the project using `run_app.bat`. Previously started localhost servers may need to be launched again after reboot.

Ubuntu and the MD backend installation are in progress and have not yet been verified as complete. No actual MD simulation has been completed. Follow-up is active again after the user completed the restart.
