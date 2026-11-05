# York (2023) — Modern Alchemical Free Energy Methods Explained

**Summary**: A rigorous perspective outlining the physical principles, statistical mechanics estimators (TI, FEP, BAR, MBAR), soft-core potentials, and diagnostic best practices of modern alchemical free energy calculations.
**Sources**: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Darrin M. York, "Modern Alchemical Free Energy Methods for Drug Discovery Explained," *ACS Physical Chemistry Au* 3(6):478–491 (2023), doi:10.1021/acsphyschemau.3c00033 (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

This perspective delivers a comprehensive, foundational explanation of state-of-the-art Alchemical Free Energy (AFE) methods. It addresses the underlying statistical mechanical estimators, transformation topologies, soft-core potential formulations, error diagnostics, and practical recommendations for utilizing AFE in lead optimization (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

## Thermodynamic cycles in alchemical design

AFE exploits the state function property of free energy by replacing the physically inaccessible process of drug binding/unbinding with computational "alchemical" transformations along a coupling parameter $\lambda \in [0, 1]$ (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).
The relative binding free energy ($\Delta\Delta G_{\text{bind}}$) between ligand A and ligand B is derived from a thermodynamic cycle:
$$\Delta\Delta G_{\text{bind}} = \Delta G_{\text{bind}}(B) - \Delta G_{\text{bind}}(A) = \Delta G_{\text{complex}}^{\text{alchem}}(A \to B) - \Delta G_{\text{water}}^{\text{alchem}}(A \to B)$$
By computing the free energy to mutate ligand A into ligand B in the bound complex and in bulk water, the net affinity difference is determined with high accuracy while bypassing real-space association kinetics (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

## Equilibrium and non-equilibrium estimators

The author details the four standard statistical mechanical formulations for extracting free energies from $\lambda$-dependent Hamiltonians $H(\lambda) = (1-\lambda) H_A + \lambda H_B$ (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf):

### 1. Thermodynamic integration (TI)
TI calculates free energy by integrating the ensemble-averaged derivative of the potential energy with respect to $\lambda$:
$$\Delta A = \int_{0}^{1} \left\langle \frac{\partial U(\lambda)}{\partial \lambda} \right\rangle_{\lambda} d\lambda$$
Evaluated numerically using quadrature or Simpson's rule across discrete intermediate $\lambda$-windows (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

### 2. Free energy perturbation (FEP / Zwanzig equation)
FEP expresses free energy through the exponential work average between states:
$$\Delta A_{i \to j} = -k_B T \ln \left\langle \exp\left(-\frac{U_j - U_i}{k_B T}\right) \right\rangle_i$$
FEP requires substantial phase space overlap between states $i$ and $j$; otherwise, the exponential average is dominated by rare, poorly sampled low-energy configurations in the target state (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

### 3. Bennett acceptance ratio (BAR) and Multistate BAR (MBAR)
- **BAR (Bennett, 1976)**: An optimal maximum-likelihood estimator that combines forward ($i \to j$) and reverse ($j \to i$) simulations to yield the minimum asymptotic variance estimator for a pair of adjacent states.
- **MBAR (Shirts & Chodera, 2008)**: A generalized extension of BAR that incorporates configuration snapshots from all sampled intermediate $\lambda$-states simultaneously. MBAR utilizes all simulation data self-consistently, eliminates histogram binning, and provides rigorous, analytical statistical error estimates (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

### 4. Non-equilibrium work (NEW) methods
NEW simulations employ rapid, non-equilibrium switching trajectories between states, recovering equilibrium free energy differences via Jarzynski's equality or Crooks Fluctuation Theorem. Combined with bidirectional maximum likelihood estimators, NEW can provide robust free energies when equilibrium sampling is difficult to maintain (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

## Transformation topologies and soft-core potentials

### 1. The endpoint catastrophe and soft-core potentials
Standard Lennard-Jones ($r^{-12}$) and Coulomb ($r^{-1}$) potentials diverge to infinity as interatomic distance $r \to 0$. When atoms are annihilated or created near $\lambda = 0$ or $\lambda = 1$, solvent or protein atoms can overlap with disappearing centers, causing numerical instability (the "endpoint singularity"). Soft-core potentials modify the radial distance by adding a $\lambda$-dependent shifting parameter:
$$r_{\text{eff}} = \left( \alpha (1-\lambda)^m \sigma^6 + r^6 \right)^{1/6}$$
which caps forces and energies at finite values as atoms vanish (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

### 2. Staged electrostatic decoupling
Bare point charges must not exist without accompanying van der Waals repulsion spheres. Consequently, the recommended transformation schedule decouples partial atomic charges linearly to zero first, followed by the soft-core transformation of Lennard-Jones parameters, preventing devastating charge clashes in the vanishing state (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

### 3. Single versus dual topology
- **Single-topology**: Atoms sharing common chemotypes are retained; only changing functional groups undergo parameter morphing. Requires fewer non-interacting dummy atoms but can introduce complex bonded constraint definitions.
- **Dual-topology**: Both ligands coexist in the simulation simultaneously, with their non-shared moieties decoupled independently. Employs dummy atoms that interact purely intramolecularly, requiring restraints to prevent dummy groups from wandering away during intermediate states (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

## Quality control, error detection, and "book-ending"

The author establishes rigorous quality-assurance criteria for AFE calculations (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf):
- **Phase space overlap matrix**: In MBAR, adjacent $\lambda$-states must maintain adequate energy overlap (typically $>0.05-0.10$). Poor overlap demands inserting intermediate $\lambda$-windows.
- **Hysteresis and estimator cross-checking**: If free energy predictions from TI, BAR, and MBAR diverge by more than their internal statistical error bounds, sampling is incomplete and hidden conformational barriers exist.
- **Cycle closure**: In relative AFE networks, summing $\Delta\Delta G$ around closed thermodynamic loops ($\sum_{\text{cycle}} \Delta\Delta G$) must equal zero within statistical error. Deviations identify specific problematic edges.
- **Indirect "book-ending" corrections**: High-level Hamiltonians (such as QM/MM or polarizable force fields) can be evaluated at the end states using thermodynamic perturbation without running exhaustive sampling at the expensive high level (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).

## Related pages
- [[binding-free-energy]]
- [[force-fields]]
- [[molecular-dynamics-simulation]]
- [[binding-thermodynamics]]
- [[enhanced-sampling]]

