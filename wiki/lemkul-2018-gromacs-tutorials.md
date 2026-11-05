# Lemkul (2019) — Suite of GROMACS Tutorials

**Summary**: Justin Lemkul's peer-reviewed LiveCoMS article presenting a comprehensive suite of seven standardized tutorials for the GROMACS simulation package, ranging from basic protein hydration to membrane bilayers, umbrella sampling, and alchemical free energy calculations.
**Sources**: `mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf`
**Last updated**: 2026-11-02

---

## Pedagogical overview and suite scope

The article documents a curated, open-access suite of seven molecular simulation tutorials designed to guide users from foundational principles to advanced biophysical free energy simulations in GROMACS (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf). The tutorials address critical entry barriers in molecular modeling, including force field selection, boundary condition artifacts, integrator settings, thermostat/barostat pairing, and command-line execution pipelines (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).

The seven tutorials span seven distinct biophysical workflows:
1. **Tutorial 1: Lysozyme in Water** — Aqueous protein setup, OPLS-AA force field, rhombic dodecahedron solvation, neutralization, $NVT$/$NPT$ equilibration, and trajectory analysis (RMSD, RMSF, radius of gyration $R_g$) (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
2. **Tutorial 2: Membrane Protein (KALP15 in DPPC)** — Insertion of a transmembrane peptide into a pre-equilibrated dipalmitoylphosphatidylcholine (DPPC) lipid bilayer using the GROMOS96 53A6 force field, grid scaling, `water_deletor.pl`, and area per lipid ($62\text{ \AA}^2$) equilibration (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
3. **Tutorial 3: Umbrella Sampling** — Computing the Potential of Mean Force (PMF) for peptide dissociation or membrane permeation along a reaction coordinate ($z$-axis), generation of series of configuration windows via steered MD, harmonic restraint biasing, and Weighted Histogram Analysis Method (WHAM) extraction via `gmx wham` (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
4. **Tutorial 4: Biphasic Systems** — Constructing self-assembling or pre-formed non-aqueous liquid–liquid interfaces (e.g., cyclohexane/water), demonstrating mutual solubility barriers and partition thermodynamics (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
5. **Tutorial 5: Protein–Ligand Complex** — Detailed parameterization of a small-molecule ligand (2-propylphenol, JZ4) using CGenFF, position restraint indexing, merged temperature coupling groups (`Protein_JZ4`), and nonbonded interaction energy decomposition (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
6. **Tutorial 6: Free Energy of Solvation** — Alchemical transformation computing the hydration free energy of methane using thermodynamic integration (TI) and Bennett Acceptance Ratio (BAR) / MBAR, with staged decoupling ($\Delta G_{\text{coul}}$ and $\Delta G_{\text{vdw}}$) and soft-core potentials (`vdw-lambdas`, `coul-lambdas`) (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
7. **Tutorial 7: Virtual Sites** — Constructing massless virtual interaction sites (VSITE) on aromatic rings and hydrogen atoms to eliminate high-frequency stretching motions, enabling an expanded $4.0\text{ fs}$ or $5.0\text{ fs}$ integration timestep for high-throughput biomolecular simulation (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).

## Key simulation best practices and guidelines

Across the suite, Lemkul establishes foundational simulation standards:
- **Unit Cell Geometry**: Rhombic dodecahedral and triclinic cells are strongly favored over cubic cells for globular proteins because they minimize the number of explicit solvent molecules by approximately $23\%$ for an identical minimum solute-to-wall distance ($1.0\text{ nm}$), markedly accelerating throughput (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
- **Thermostat Grouping**: Small solutes or bound ligands must not be coupled to independent thermostats or grouped with bulk solvent; they must be merged with the macromolecule (`tc-grps = Protein_Ligand Water_and_ions`) to prevent spurious thermal gradients (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
- **Alchemical Staging**: In free energy calculations, electrostatic interactions must be completely uncoupled before Lennard-Jones interactions are modified to avoid the "endpoint catastrophe" where charges overlap at zero separation, or soft-core potentials must be applied with `sc-power = 1` and `sc-alpha = 0.5` (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
- **Equilibration Checks**: Trajectories must be verified for thermodynamic stability before data acquisition: constant energy in $NVE$, converged temperature around $T_0$ in $NVT$, and stable density and box dimensions in $NPT$ (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).

## Trajectory analysis and post-processing protocols

The suite highlights post-processing utilities in GROMACS:
- Periodic boundary condition (PBC) handling requires a two-step procedure: centering the solute with `gmx trjconv -center -pbc mol -ur compact`, followed by rotational and translational fitting via `-fit rot+trans` (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
- Standard biophysical trajectory analyses include `gmx rms` (structural deviation), `gmx rmsf` (residue flexibility/B-factor correlation), `gmx gyrate` (mass compactness), `gmx hbond` (hydrogen-bond persistence), and `gmx wham` (potential of mean force reconstruction) (source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).

## Related pages
- [[gromacs-protein-ligand-protocol]]
- [[lemkul-gromacs-protein-ligand-tutorial]]
- [[molecular-dynamics-simulation]]
- [[umbrella-sampling]]
- [[alchemical-free-energy-methods]]
- [[force-fields]]
- [[simulation-and-benchmark-best-practices]]
- [[hydrogen-mass-repartitioning]]

