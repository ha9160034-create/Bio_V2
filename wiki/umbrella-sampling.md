# Umbrella Sampling

**Summary**: A collective variable-based free energy method that restrains a molecular system across discrete spatial windows using harmonic bias potentials to reconstruct potential of mean force profiles and standard binding free energies.
**Sources**: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulations and drug discovery.md; Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf
**Last updated**: 2026-11-01

---

## Physical principles and harmonic bias formulation

Umbrella sampling, introduced by Torrie and Valleau (1977), is a foundational enhanced sampling method designed to compute free energy profiles along defined reaction coordinates by overcoming high kinetic barriers that prevent ergodic sampling in unbiased [[molecular-dynamics-simulation]] (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

In an unbiased simulation, the probability $P(\boldsymbol{\xi})$ of observing a system at a specific configuration coordinate $\boldsymbol{\xi}$ is proportional to the Boltzmann factor $\exp(-F(\boldsymbol{\xi})/k_B T)$ (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). High free-energy barriers ($F \gg k_B T$) are sampled with vanishingly low probability. Umbrella sampling surmounts this limitation by adding an external, artificial biasing potential—most commonly a harmonic restraint—that bridges the high-energy regions (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):
$$U_{\text{bias}, i}(\mathbf{r}) = \frac{k_i}{2} \left[ \xi(\mathbf{r}) - \xi_{0, i} \right]^2$$
where $\xi(\mathbf{r})$ is the reaction coordinate (e.g., center-of-mass distance between ligand and receptor binding site), $\xi_{0, i}$ is the center of the $i$-th umbrella window, and $k_i$ is the harmonic force constant (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Stratified sampling and window design

The reaction coordinate spanning the bound complex to the fully dissociated state is discretized into a series of $N_w$ overlapping umbrella windows (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):
1. **Starting conformations**: Initial atomic configurations for each window are generated via Steered Molecular Dynamics (SMD), pulling the ligand slowly out of the pocket along $\xi$, or through high-temperature trajectory seeding (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulations and drug discovery.md).
2. **Equilibration and production**: Each window is simulated independently with its harmonic bias restraint active, sampling local configurations around $\xi_{0, i}$ (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
3. **Histogram overlap criterion**: The spacing between window centers ($\Delta \xi_0 = \xi_{0, i+1} - \xi_{0, i}$) and the force constants ($k_i$) must be calibrated so that adjacent configuration histograms exhibit substantial statistical overlap (typically $>15\text{--}20\%$ overlap area) (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). If histograms do not overlap, unbiasing algorithms fail to connect adjacent free energy windows, producing catastrophic numerical offsets (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Reconstructing the potential of mean force (PMF)

The biased probability distribution $P_i^b(\xi)$ collected from the $i$-th window is related to the unbiased probability distribution $P^u(\xi)$ by:
$$P^u(\xi) = P_i^b(\xi) \exp\left( \frac{U_{\text{bias}, i}(\xi)}{k_B T} \right) \left\langle \exp\left( -\frac{U_{\text{bias}, i}(\xi)}{k_B T} \right) \right\rangle$$
(source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). Two primary statistical estimators are utilized to combine multiple overlapping windows:

### 1. Weighted histogram analysis method (WHAM)
WHAM calculates the optimal unbiased distribution by solving a set of self-consistent non-linear equations that minimize the statistical variance across all windows:
$$P^u(\xi) = \frac{\sum_{i=1}^{N_w} N_i P_i^b(\xi)}{\sum_{i=1}^{N_w} N_i \exp\left( -\frac{U_{\text{bias}, i}(\xi) - F_i}{k_B T} \right)}$$
$$F_i = -k_B T \ln \left( \int P^u(\xi) \exp\left( -\frac{U_{\text{bias}, i}(\xi)}{k_B T} \right) d\xi \right)$$
where $N_i$ is the number of trajectory frames in window $i$, and $F_i$ is the relative free energy offset of window $i$ (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). The Potential of Mean Force is then obtained directly as $W(\xi) = -k_B T \ln P^u(\xi)$ (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### 2. Multistate Bennett acceptance ratio (MBAR)
MBAR provides a modern, unbinned maximum-likelihood estimator that eliminates histogram binning artifacts (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). By evaluating configuration energies across all sampled states simultaneously, MBAR extracts the lowest-variance free energy curve and computes analytical asymptotic confidence intervals (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

## Converting PMF to standard binding free energy ($\Delta G^\circ$)

A one-dimensional PMF profile $W(\xi)$ does not equal the standard binding free energy ($\Delta G^\circ$) directly; it must be corrected for standard-state volume and rotational/translational restraints (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):
1. **Partition function integration**: The PMF must be integrated over the bound state region ($\text{site}$) and the bulk solvent region ($\text{bulk}$):
$$K_a = \frac{\int_{\text{site}} \exp\left( -\frac{W(\xi)}{k_B T} \right) d\xi}{\int_{\text{bulk}} \exp\left( -\frac{W(\xi)}{k_B T} \right) d\xi}$$
(source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
2. **Standard-state volume correction**: The standard state is defined at an effective reactant concentration of $C^\circ = 1\text{ M}$, corresponding to a standard reference volume of:
$$V_0 = \frac{1}{C^\circ N_A} = 1661\text{ \AA}^3$$
(source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). Moving from the volume actually sampled in the simulation box ($V_{\text{box}}$) to $V_0$ introduces a mandatory analytical correction $\Delta G_{\text{vol}} = -k_B T \ln(V_{\text{box}} / V_0)$ (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
3. **Restraint contributions**: When harmonic restraints are applied orthogonal to the reaction coordinate (e.g., cylindrical or Boresch orientation restraints to keep the ligand aligned during unbinding), the free energy cost of removing these restraints in the bulk state ($\Delta G_{\text{restraint}}$) must be added (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). Without these corrections, PMF depths deviate from experimental calorimetry by several kcal/mol (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Comparison with other enhanced sampling methods

| Feature | Umbrella Sampling | Metadynamics | Steered MD (SMD) |
|---|---|---|---|
| **Sampling mode** | Equilibrium in stratified windows | Non-equilibrium history-dependent bias | Non-equilibrium forced pulling |
| **Free energy estimator** | WHAM / MBAR | Negative of accumulated bias ($-V$) | Jarzynski equality / Crooks theorem |
| **Convergence control** | High (controlled window histograms) | High with Well-Tempered (WT-MetaD) | Poor if pulling speed is rapid |
| **Computational cost** | High (requires 20–60 separate windows) | Medium-High (continuous trajectory) | Low-Medium (replicated pulling runs) |
| **Kinetic parameters** | Indirect (requires diffusive Kramers modeling) | Direct via Infrequent Metadynamics (InMetaD) | Qualitative rupture force profiling |

Sources: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf.

## Related pages
- [[enhanced-sampling]]
- [[binding-free-energy]]
- [[metadynamics]]
- [[molecular-dynamics-simulation]]
- [[alchemical-free-energy-methods]]
- [[binding-thermodynamics]]
