# Genheden & Ryde (2015) — The MM/PBSA and MM/GBSA Methods

**Summary**: An authoritative review detailing the theoretical foundation, parameter calibration, precision limits, and practical best practices of the MM/PBSA and MM/GBSA endpoint binding free energy methods.
**Sources**: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Samuel Genheden and Ulf Ryde, "The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities," *Expert Opinion on Drug Discovery* 10(5):449–461 (2015), doi:10.1517/17460441.2015.1032936 (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

The review evaluates the Molecular Mechanics/Poisson–Boltzmann Surface Area (MM/PBSA) and Molecular Mechanics/Generalized Born Surface Area (MM/GBSA) methodologies. These endpoint approaches occupy an intermediate niche in computer-aided drug design: more physically grounded and accurate than empirical docking scoring functions, but orders of magnitude faster than strict alchemical free energy perturbation (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

## Thermodynamic formulation

The binding free energy is approximated using a thermodynamic cycle:
$$\Delta G_{\text{bind}} = \langle G_{\text{complex}} \rangle - \langle G_{\text{receptor}} \rangle - \langle G_{\text{ligand}} \rangle$$
where the free energy of each state is partitioned into gas-phase molecular mechanics, continuum solvation, and conformational entropy:
$$G = E_{\text{MM}} + G_{\text{solv}} - T S_{\text{conf}}$$
$$E_{\text{MM}} = E_{\text{bonded}} + E_{\text{vdw}} + E_{\text{elec}}$$
$$G_{\text{solv}} = G_{\text{polar}} + G_{\text{non-polar}}$$
(source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

## Single-trajectory versus three-trajectory approach

A fundamental distinction highlighted by the authors is whether ensembles are sampled from one or three separate MD simulations (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf):
- **Single-trajectory approach (1-trajectory)**: Trajectories are generated solely for the protein–ligand complex. Snapshots of the free receptor and free ligand are derived by stripping coordinates from the complex trajectory. Crucially, the internal bonded terms ($E_{\text{bond}}$, $E_{\text{angle}}$, $E_{\text{dihedral}}$) cancel exactly ($\Delta E_{\text{bonded}} = 0$), drastically reducing statistical variance and noise.
- **Three-trajectory approach (3-trajectory)**: Complex, apo-receptor, and free ligand are simulated independently. Although theoretically capturing receptor and ligand reorganization energies, large fluctuations in internal force-field energies yield standard deviations on the order of 10–20 kcal/mol, completely masking affinity differences unless impractical sampling lengths are used.
- **Conclusion**: The 1-trajectory approach is universally recommended for virtual screening and congeneric lead ranking (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

## Solvation energy models

### 1. Polar solvation ($G_{\text{polar}}$)
- **Poisson–Boltzmann (PB)**: Solves the linear or non-linear PB differential equation numerically on a cubic grid. PB provides rigorous electrostatic continuum modeling but is computationally heavier.
- **Generalized Born (GB)**: An analytical approximation of PB based on effective Born radii of atoms. GB is computationally fast, differentiable, and generally performs comparably to PB in relative ranking tasks (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

### 2. Internal solute dielectric constant ($\epsilon_{\text{int}}$)
The authors emphasize that selecting $\epsilon_{\text{int}}$ is critical (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf):
- $\epsilon_{\text{int}} = 1$ corresponds to a rigid gas-phase vacuum solute; it tends to heavily exaggerate electrostatic interactions.
- $\epsilon_{\text{int}} = 2$ or $4$ models electronic polarization and minor side-chain relaxation. Ryde's systematic benchmarks show that $\epsilon_{\text{int}} = 2$ or $4$ often yields higher correlations with experimental affinities than $\epsilon_{\text{int}} = 1$, especially for charged ligands, solvent-exposed active sites, or metalloenzymes (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

### 3. Non-polar solvation ($G_{\text{non-polar}}$)
- **SASA-only model**: $\Delta G_{\text{non-polar}} = \gamma \cdot \text{SASA} + \beta$, where $\gamma$ is surface tension and $\beta$ an offset.
- **Cavitation + dispersion model**: Separates non-polar free energy into unfavorable cavity formation ($G_{\text{cav}}$) and favorable solute–solvent van der Waals attractive dispersion ($G_{\text{disp}}$), which provides a more realistic physical representation of hydrophobic burial (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

## The entropy conundrum ($-T\Delta S$)

The loss of conformational entropy upon binding is typically computed via Normal Mode Analysis (NMA) or Quasi-Harmonic Analysis (QHA) (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf):
- NMA requires extensive conjugate-gradient energy minimization of snapshots in vacuum, which can distort the geometry and introduce substantial numerical noise ($\sim 2-5\text{ kcal/mol}$).
- **Key benchmark insight**: In many drug discovery benchmarks evaluating relative affinities of congeneric series, **omitting the $-T\Delta S$ term entirely improves experimental correlation**. Because entropic variation among closely related analogs is small, the computational uncertainty of NMA exceeds the real signal, degrading predictive power (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

## Best practices for precision and simulation length

- **Multiple short trajectories**: Running 5 to 10 independent short MD simulations (e.g., 1–2 ns each) with randomized initial velocities or varied initial water/protonation states provides greater precision and lower standard error than a single long trajectory (e.g., 20–50 ns), which frequently suffers from irreversible force-field drift into non-native conformations (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).
- **Snapshot frequency**: Sampling 100 to 500 uncorrelated snapshots is typically sufficient for convergence in the 1-trajectory framework (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

## Related pages
- [[binding-free-energy]]
- [[molecular-dynamics-simulation]]
- [[virtual-screening]]
- [[force-fields]]
- [[water-and-solvation]]
- [[scoring-functions]]

