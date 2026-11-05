# Lazim, Suh & Choi (2020) — Advances in MD and Enhanced Sampling

**Summary**: A review of modern advances in molecular dynamics, classifying enhanced sampling into CV-based and non-CV-based methods, and evaluating their application to GPCRs, kinases, and protein complexes.
**Sources**: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Raudah Lazim, Donghyuk Suh, and Sun Choi, "Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems," *International Journal of Molecular Sciences* 21(17):6339 (2020), doi:10.3390/ijms21176339 (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

The review evaluates computational advancements in atomistic molecular dynamics (MD) simulations, detailing algorithmic optimizations such as hydrogen mass repartitioning, systematically classifying enhanced sampling methods into collective variable (CV)-dependent and non-CV-dependent architectures, and surveying applications in computer-aided drug discovery (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

## Technical advancements in classical simulation

### Hydrogen mass repartitioning (HMR)
The authors detail Hydrogen Mass Repartitioning as a high-impact efficiency enhancement for routine biomolecular MD (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf):
- **Mechanism**: The mass of hydrogen atoms is increased from ~1.008 amu to ~3.0 amu by transferring mass from the bonded heavy atom (carbon, nitrogen, or oxygen), preserving the overall molecular mass.
- **Biophysical consequence**: Increasing the mass of hydrogen slows down the highest-frequency bond stretching vibrations (e.g., C–H, N–H, O–H bonds), allowing the numerical integration timestep to be increased from 2.0 fs to **4.0 fs**.
- **Thermodynamic validity**: Benchmark evaluations on hen egg white lysozyme (HEWL) and phospholipid bilayers confirm that HMR grants a **two-fold increase in simulation throughput** without inducing statistically significant deviations in thermodynamic free energy profiles, membrane area-per-lipid, or structural order parameters (Hopkins et al., 2015; Balusek et al., 2019) (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

## Classification of enhanced sampling methodologies

The authors categorize enhanced sampling strategies into two overarching paradigms (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf):

### 1. Collective variable (CV)-based methods
These methods require the explicit selection of geometric order parameters (e.g., center-of-mass distances, dihedral angles, coordination numbers, or RMSD):
- **Umbrella Sampling (US)**: Restrains the system along a reaction coordinate using harmonic umbrella potentials; the full Potential of Mean Force (PMF) is unbiased using the Weighted Histogram Analysis Method (WHAM) or Dynamic WHAM (D-WHAM).
- **Metadynamics variants**:
  - *Standard Metadynamics*: Deposits repulsive Gaussians along selected CVs.
  - *Well-Tempered Metadynamics (WT-MetaD)*: Reduces Gaussian hill height progressively as a function of deposited potential, ensuring smooth convergence to the free energy minimum.
  - *Bias-Exchange Metadynamics (BE-MetaD)*: Simulates multiple replicas simultaneously, where each replica biases a distinct CV, exchanging configurations via a Metropolis criterion to sample complex multidimensional landscapes.
- **Adaptive Biasing Force (ABF)**: Directly estimates the mean thermodynamic force acting along the CV and applies an equal, opposite cancellation force to induce a uniform diffusive walk without free-energy barriers.
- **Steered MD (SMD)**: Applies constant velocity or constant force springs to pull ligands out of pockets or unfold domains, reconstructing work profiles via Jarzynski's equality (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

### 2. Non-CV-based methods
These methods accelerate state transitions across the entire potential energy surface without requiring prior knowledge or manual selection of reaction coordinates:
- **Accelerated MD (aMD)**: Adds a non-negative boost potential to the total or torsional potential when the system potential drops below a reference energy $E$.
- **Gaussian Accelerated MD (GaMD)**: Miao et al. refined aMD by imposing harmonic boost potentials calibrated such that the boosted potential energy follows a near-Gaussian distribution. This crucial innovation enables accurate free-energy unbiasing via cumulant expansion to the second order, resolving the severe reweighting noise that plagued original aMD.
- **Replica Exchange (T-REMD / H-REMD)**: Exchanges trajectories across temperature or Hamiltonian spaces to cross orthogonal barriers without geometric constraints (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

## Applications to drug discovery target classes

The authors review how these advanced techniques resolve challenging pharmacotherapeutic targets (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf):
- **G Protein-Coupled Receptors (GPCRs)**: Capturing active, intermediate, and inactive states (e.g., ionic lock formation between DRY motif $Arg^{3.50}$ and $Glu^{6.30}$), allosteric modulator binding in lipid-facing crevices, and G-protein coupling.
- **Kinases**: Investigating DFG-in / DFG-out activations, $\alpha$C-helix conformations, and activation loop phosphorylation states using GaMD and WT-MetaD.
- **Protein–Protein & Protein–Peptide Interfaces**: Combining rigid-body docking (ClusPro PeptiDock) with flexible MD refinement (Rosetta FlexPepDock) and MM-PBSA/GBSA to characterize dynamic interface hot spots (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

## Related pages
- [[enhanced-sampling]]
- [[molecular-dynamics-simulation]]
- [[binding-thermodynamics]]
- [[force-fields]]
- [[binding-kinetics]]
- [[protein-ligand-interactions]]

