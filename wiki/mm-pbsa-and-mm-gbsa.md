# MM-PBSA and MM-GBSA Methods

**Summary**: Endpoint binding free energy methods that combine molecular mechanics force fields with continuum Poisson–Boltzmann or Generalized Born solvation to rescore docking poses and estimate lead affinity.
**Sources**: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf; Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md
**Last updated**: 2026-11-01

---

## Methodological position in CADD

Molecular Mechanics/Poisson–Boltzmann Surface Area (MM-PBSA) and Molecular Mechanics/Generalized Born Surface Area (MM-GBSA) occupy an intermediate tier in computer-aided drug discovery:
- Significantly more physically rigorous and accurate than fast empirical or knowledge-based docking [[scoring-functions]] because they incorporate dynamic conformational sampling, molecular mechanics force fields, and continuum solvent electrostatics (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- Orders of magnitude faster than rigorous [[alchemical-free-energy-methods]] (FEP/TI), making them computationally tractable for rescoring hundreds to thousands of molecules in [[virtual-screening]] (source: Application of MM-PBSA Methods in Virtual Screening.pdf; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).
- However, as reviewed by Bruce et al. (2018), MM-PB(GB)SA is an approximate endpoint method lacking systematic mathematical convergence; longer trajectories do not necessarily yield higher correlation, and absolute affinities carry substantial errors ($\sim 5\text{ kcal/mol}$), making the method suited primarily for relative ranking and hit rescoring (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). Trajectory generation can be substantially accelerated using [[hydrogen-mass-repartitioning]] (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Thermodynamic formulation

The binding free energy is defined by the standard thermodynamic cycle:
$$\Delta G_{\text{bind}} = \langle G_{\text{complex}} \rangle - \langle G_{\text{receptor}} \rangle - \langle G_{\text{ligand}} \rangle$$
where each free energy state is calculated as:
$$G = E_{\text{MM}} + G_{\text{solv}} - T S_{\text{conf}}$$
$$E_{\text{MM}} = E_{\text{bonded}} + E_{\text{vdw}} + E_{\text{elec}}$$
$$G_{\text{solv}} = G_{\text{polar}} + G_{\text{non-polar}}$$
(source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

## Single-trajectory versus three-trajectory approach

A defining methodological choice in MM-PB(GB)SA calculations is trajectory sampling (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf):
- **Single-trajectory approach (1-trajectory)**: Only the protein–ligand complex is simulated via [[molecular-dynamics-simulation]]. Coordinate snapshots for the free receptor and free ligand are derived by stripping atoms from the complex trajectory. Because the conformations are identical, internal bonded terms ($E_{\text{bond}}$, $E_{\text{angle}}$, $E_{\text{dihedral}}$) cancel exactly ($\Delta E_{\text{bonded}} = 0$). This eliminates massive thermal noise and dramatically lowers standard errors (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).
- **Three-trajectory approach (3-trajectory)**: Simulates the complex, unbound receptor, and unbound ligand independently. While theoretically accounting for receptor and ligand conformational reorganization upon binding, the fluctuations in bonded terms produce enormous standard deviations (typically 10–20 kcal/mol), completely destroying correlation with experimental affinity (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).
- **Consensus**: The single-trajectory approach is almost universally recommended for virtual screening rescoring and congeneric lead ranking (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf).

## Continuum solvation models

### 1. Polar solvation ($G_{\text{polar}}$)
- **Poisson–Boltzmann (PB)**: Numerically solves the PB equation on a spatial grid, modeling the solvent as a high-dielectric continuum ($\epsilon \approx 80$) and the solute as a low-dielectric cavity ($\epsilon_{\text{int}}$).
- **Generalized Born (GB)**: An analytical approximation of PB based on atom pairwise distances and effective Born radii. GB is significantly faster than PB, analytical, and well suited for rescoring large compound batches in [[virtual-screening]] (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf).

### 2. Internal solute dielectric constant ($\epsilon_{\text{int}}$)
Selecting $\epsilon_{\text{int}}$ is crucial for balancing electrostatic interaction and desolvation penalties (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf):
- $\epsilon_{\text{int}} = 1$: Rigid vacuum representation; often overestimates electrostatic interactions and salt-bridge strength.
- $\epsilon_{\text{int}} = 2$ or $4$: Accounts for electronic polarization and localized side-chain relaxation. Benchmarks indicate that $\epsilon_{\text{int}} = 2-4$ frequently yields higher correlation with experimental binding affinities, particularly for charged ligands, solvent-exposed active sites, or metalloproteins (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf).

### 3. Non-polar solvation ($G_{\text{non-polar}}$)
- **SASA-only model**: $\Delta G_{\text{non-polar}} = \gamma \cdot \text{SASA} + \beta$, scaling linearly with solvent accessible surface area.
- **Cavitation + dispersion model**: Explicitly divides non-polar free energy into unfavorable cavity formation ($G_{\text{cav}}$) and favorable solute–solvent van der Waals attractive dispersion ($G_{\text{disp}}$), providing a more realistic physical representation of [[water-and-solvation]] (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf).

## The entropy dilemma ($-T\Delta S$)

Conformational entropy loss upon binding is traditionally estimated using Normal Mode Analysis (NMA) or Quasi-Harmonic Analysis (QHA) (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf):
- NMA requires aggressive in vacuo energy minimization of snapshots, which can introduce severe structural distortion and large numerical errors ($\sim 2-5\text{ kcal/mol}$).
- **The omission rule**: In congeneric series ranking and virtual screening, **omitting the $-T\Delta S$ term entirely often improves experimental correlation**. Because entropic variations among structurally related ligands are small, the computational noise of NMA exceeds the true physical signal (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf).
- **Interaction entropy (IE)**: Duan et al. formulated interaction entropy to extract entropy directly from gas-phase interaction energy fluctuations. While effective in some targets (such as Bcl-2 family members), it can introduce noise in flexible active sites (such as EGFR), demonstrating that entropy calculation remains target-dependent (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

## Deployment in virtual screening and lead optimization

In modern hierarchical [[virtual-screening]] workflows, MM-PBSA and MM-GBSA serve as an intermediate post-docking filter (source: Application of MM-PBSA Methods in Virtual Screening.pdf; Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf):
1. **Primary docking**: Filters libraries from $10^6$ compounds to the top 500–2,000 poses.
2. **MD-based MM-PBSA rescoring**: 0.5 to 5 ns explicit solvent MD simulation relaxes the complex, filtering out false positives that detach or suffer heavy desolvation penalties, boosting 1% Enrichment Factors (EF1%) from ~1.2 to >5.0 across targets like COX-2, Thrombin, and BACE1 (source: Application of MM-PBSA Methods in Virtual Screening.pdf).
3. **Hit selection**: Top 50–100 rescored candidates proceed to visual inspection, structural clustering, and experimental validation (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

## Related pages
- [[binding-free-energy]]
- [[virtual-screening]]
- [[molecular-dynamics-simulation]]
- [[hydrogen-mass-repartitioning]]
- [[alchemical-free-energy-methods]]
- [[scoring-functions]]
- [[water-and-solvation]]
- [[force-fields]]
- [[benchmark-drug-targets]]
- [[genheden-2015-mm-pbsa-gbsa]]
- [[poli-2020-mm-pbsa-virtual-screening]]
