## Code & Pipeline Development Rules (Wiki as Mandatory Brain)

When developing, extending, or maintaining the Python code for Molecular Docking and Molecular Dynamics (MD):

1. **Mandatory Wiki Lookup Before Coding**:
   - NEVER invent simulation parameters, CLI commands, or physical workflows from general training memory.
   - Before writing or modifying any script in the `md/` or `docking/` packages, read the relevant concept and protocol pages in `wiki/`:
     - Docking & Pose Selection: `wiki/autodock-vina.md`, `wiki/docking-failure-modes.md`, `wiki/deep-learning-docking.md`.
     - MD Engine & Setup: `wiki/gromacs-protein-ligand-protocol.md`, `wiki/molecular-dynamics-simulation.md`, `wiki/simulation-and-benchmark-best-practices.md`.
     - Small-Molecule Topology: `wiki/open-force-field-and-ligand-parameterization.md`, `wiki/force-fields.md`.
     - Post-MD Analysis & Rescoring: `wiki/mm-pbsa-and-mm-gbsa.md`, `wiki/binding-free-energy.md`.

2. **Strict Biophysical Implementation Guardrails**:
   - Unit Cell: Default to rhombic dodecahedral geometry (`-bt dodecahedron -d 1.0`) to save ~23% solvent volume.
   - Restraints & Thermostats: Restrain ligand heavy atoms during solvent equilibration, and merge protein and ligand into a single thermostat coupling group (`tc-grps = Protein_Ligand Water_and_ions`).
   - Ensemble Validity: Use Bussi velocity-rescaling (`v-rescale`) or Nosé–Hoover for temperature coupling. NEVER use the Berendsen thermostat/barostat in production simulations.
   - Trajectory Processing: Always apply periodic boundary condition (PBC) re-wrapping (`trjconv -pbc mol -ur compact`) before calculating RMSD, RMSF, or H-bonds.
   - Traceability: Add a top docstring/comment in every generated Python module citing the specific `wiki/*.md` pages used to design its workflow.