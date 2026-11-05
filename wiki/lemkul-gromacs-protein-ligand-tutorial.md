# Lemkul (Tutorial) — GROMACS Protein–Ligand Complex

**Summary**: Justin Lemkul's authoritative tutorial detailing the complete end-to-end simulation workflow for a protein–ligand complex (T4 lysozyme L99A/M102Q bound to 2-propylphenol, JZ4) in GROMACS using the CHARMM36 all-atom force field and CGenFF.
**Sources**: `gromacs_protein_ligand_complex_tutorial_lemkul.md`
**Last updated**: 2026-11-02

---

## Overview and model system

The tutorial details the setup, simulation, and trajectory analysis of a small-molecule ligand bound to a globular protein in [[gromacs-protein-ligand-protocol]] (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). The model biological system comprises the mutant bacteriophage **T4 lysozyme L99A/M102Q** (PDB code `3HTB`, solved at high crystallographic resolution) complexed with the neutral synthetic ligand **2-propylphenol** (PDB chemical identifier `JZ4`, 22 atoms) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

Because small-molecule ligands are absent from standard GROMACS residue topology libraries (`.rtp`), the system must be decoupled into separate protein and ligand preparation pipelines prior to final recombination in the master system topology (`topol.top`) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Topology generation and ligand parameterization

The protein coordinates are cleaned of crystallographic solvent and crystallization agents (`PO4` phosphate and `BME` $\beta$-mercaptoethanol) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). Passing `3HTB_clean.pdb` to `gmx pdb2gmx` generates the protein coordinate file (`3HTB_processed.gro`), position restraint file (`posre.itp`), and initial topology using the **CHARMM36 all-atom force field** and CHARMM-modified **TIP3P** water model with neutral terminal patches (`NH3+` and `COO-`), yielding a net protein charge of **+6e** (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

Small-molecule parameterization follows the CHARMM General Force Field (**CGenFF**) protocol (source: gromacs_protein_ligand_complex_tutorial_lemkul.md):
- Explicit hydrogen atoms are added via Avogadro to satisfy valence requirements, exported as a Sybyl `.mol2` file, and formatted with unified residue IDs (`1 JZ4`) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Bond ordering is strictly standardized using `sort_mol2_bonds.pl` prior to submission to the CGenFF server (`cgenff.com` or `paramchem.org`) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Parametrization reliability is governed by CGenFF penalty scores: penalties $< 10$ denote highly reliable parameters ready for production MD (for JZ4, penalties are $< 1$); scores between $10\text{ and }50$ warrant basic QM validation (e.g., torsional potential energy surface scans); and scores $> 50$ necessitate extensive manual QM re-parameterization (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- The converted GROMACS ligand topology (`jz4_gmx.itp`) is inserted into `topol.top` immediately after `#include "charmm36.ff/forcefield.itp"` and before the protein `[ moleculetype ]` block to ensure custom dihedral parameters are registered before molecule construction (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Solvation, electrostatics, and force-switch parameters

The complex is placed inside a **rhombic dodecahedral unit cell** with a minimum solute-to-box distance of $1.0\text{ nm}$ (`-bt dodecahedron -d 1.0`), which saves approximately $23\%$ in water volume compared to a simple cubic box of equivalent periodic buffer distance (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). The cell is solvated (`gmx solvate`), and 6 chloride counterions (`CL`) are introduced via `gmx genion` to replace solvent water molecules and neutralize the net $+6e$ charge (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

Energy minimization is performed using steepest descents until $F_{\text{max}} < 1000\text{ kJ/(mol}\cdot\text{nm)}$ (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). Under CHARMM36, Lennard-Jones van der Waals interactions require a continuous force-switching function (`vdw-modifier = force-switch`) between `rvdw-switch = 1.0` nm and `rvdw = 1.2` nm, coupled to Particle Mesh Ewald (`coulombtype = PME`) with `rcoulomb = 1.2` nm (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Restraints, temperature coupling, and production MD

During solvent relaxation, both the protein and the ligand heavy atoms are restrained to prevent destabilization of the crystallographic binding pose (source: gromacs_protein_ligand_complex_tutorial_lemkul.md):
- Heavy-atom restraints for the ligand (`posre_jz4.itp`, force constant $1000\text{ kJ/(mol}\cdot\text{nm}^2)$) are generated via `gmx make_ndx` and `gmx genrestr` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Merged thermostat groups (`tc-grps = Protein_JZ4 Water_and_ions`) are created in `index.ndx` using `1 | 13` to avoid artificial kinetic temperature oscillations caused by coupling a 22-atom solute to an independent thermostat or lumping it into bulk solvent (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- The system undergoes two sequential 100-ps equilibration phases: canonical ($NVT$, $300\text{ K}$) using modified Berendsen/v-rescale thermostats, followed by isothermal–isobaric ($NPT$, $1\text{ bar}$) with the Berendsen or Parrinello–Rahman barostat (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- A 10-ns unrestrained production MD run is executed with a $2.0\text{ fs}$ timestep using leap-frog integration and LINCS holonomic bond constraints (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Trajectory analysis and interaction metrics

Post-simulation trajectory processing resolves periodic boundary condition (PBC) artifacts:
- Molecules broken across periodic cell edges are re-wrapped into a whole, compact representation centered on the protein solute via `gmx trjconv -pbc mol -ur compact -center` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Rotational and translational least-squares fitting against the protein `Backbone` is executed via `gmx trjconv -fit rot+trans` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Hydrogen bonding between the phenolic hydroxyl of JZ4 (`OAB`/`H12`) and the sidechain carbonyl oxygen (`OE1`) of **Gln102** is quantified: the average donor–acceptor distance is $0.31 \pm 0.05\text{ nm}$ (satisfying the geometric cutoff $\le 0.35\text{ nm}$), and the deviation angle is $23.3 \pm 9.4^\circ$ (satisfying the linearity cutoff $\le 30^\circ$) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Ligand binding stability is measured by heavy-atom RMSD fitted to the protein backbone, which plateaus at $\sim 0.1\text{ nm } (1.0\text{ \AA})$ (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Nonbonded short-range interaction energy decomposition via `gmx mdrun -rerun` yields $\text{Coul-SR} = -20.5 \pm 7.4\text{ kJ/mol}$ and $\text{LJ-SR} = -99.1 \pm 7.2\text{ kJ/mol}$ (total $-119.6 \pm 10.3\text{ kJ/mol}$), illustrating gas-phase electrostatic and van der Waals contacts while noting that raw interaction energies omit desolvation and entropic penalties required for true free energies of binding (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Related pages
- [[gromacs-protein-ligand-protocol]]
- [[molecular-dynamics-simulation]]
- [[force-fields]]
- [[open-force-field-and-ligand-parameterization]]
- [[simulation-and-benchmark-best-practices]]
- [[ligand-preparation-and-protonation]]
- [[receptor-preparation]]
- [[lemkul-2018-gromacs-tutorials]]

