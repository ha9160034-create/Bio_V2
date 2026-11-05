# Ahmad et al. (2022) — Enhanced Sampling for Ligand Binding Kinetics

**Summary**: A review evaluating enhanced-sampling molecular dynamics methodologies for predicting drug unbinding rate constants ($k_{\text{off}}$) and residence times.
**Sources**: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Katya Ahmad, Andrea Rizzi, Riccardo Capelli, Davide Mandelli, Wenping Lyu, and Paolo Carloni, "Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics: Current Status and Perspective," *Frontiers in Molecular Biosciences* 9:899805 (2022), doi:10.3389/fmolb.2022.899805 (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

The review examines atomistic enhanced-sampling molecular dynamics (MD) methodologies designed to compute binding kinetics—specifically the association rate constant ($k_{\text{on}}$), dissociation rate constant ($k_{\text{off}}$), and drug residence time ($\tau = 1/k_{\text{off}}$). Because drug unbinding takes place over timescales ranging from seconds to hours ($k_{\text{off}} \approx 10^{-4}-10^2\text{ s}^{-1}$), classical unbiased simulations cannot capture spontaneous dissociation without specialized acceleration (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

## Pharmacological significance of binding kinetics

The authors emphasize that static equilibrium affinity ($\Delta G_{\text{bind}}$ or $K_d = k_{\text{off}}/k_{\text{on}}$) is often an incomplete predictor of in vivo pharmacological efficacy (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf):
- Under open, non-equilibrium physiological conditions, drug concentration fluctuates due to pharmacokinetic clearance.
- **Drug residence time ($\tau = 1/k_{\text{off}}$)** dictates the duration of target engagement: ligands with slow off-rates sustain therapeutic blockade long after systemic drug levels have cleared, conferring prolonged efficacy, target selectivity, and reduced toxicity (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

## Enhanced sampling methodologies for kinetic rate estimation

The authors systematically review the principal computational frameworks for estimating $k_{\text{off}}$ (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf):

### 1. Infrequent metadynamics (InMetaD)
Developed by Tiwary and Parrinello (2013), InMetaD adapts standard metadynamics to extract unperturbed kinetic rates (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf):
- **Deposition pacing**: Gaussian bias hills are deposited at a frequency significantly slower than the transition-state crossing time ($\tau_{\text{dep}} \gg \tau_{\text{trans}}$).
- **Physical principle**: By allowing the system to thermalize between bias additions, the repulsive potential fills only the bound free energy well without depositing bias on the transition state barrier.
- **Poisson statistics**: The unbinding transition times across independent replicas must obey a single-exponential Poisson distribution, verified via the Kolmogorov–Smirnov test ($p > 0.05$).
- **Acceleration factor**: The true unbinding time is recovered by scaling simulation time $t_{\text{sim}}$ by the time-averaged acceleration factor:
  $$\alpha = \left\langle \exp\left(\frac{V(s, t)}{k_B T}\right) \right\rangle$$
  (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

### 2. Frequency-adaptive metadynamics and OPES-flooding
- **Insubria frequency adaptation**: Dynamically scales the bias deposition rate depending on the position in CV space, accelerating sampling in flat regions while protecting barrier crossing integrity.
- **OPES-flooding**: An adaptation of On-the-fly Probability Enhanced Sampling (OPES) that fills the bound basin rapidly and stably, providing accelerated transition time estimation with reduced parameter tuning (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

### 3. $\tau$-Random acceleration molecular dynamics ($\tau$-RAMD)
Developed by Kokh, Wade, and co-workers (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf):
- **Mechanism**: A random force vector is applied to the ligand center of mass. If the ligand fails to move past a threshold distance within a short evaluation window, the force direction is randomly reassigned.
- **Application**: Generates complete dissociation trajectories in nanoseconds. When averaged over 50–100 runs, computed relative unbinding times correlate strongly with experimental residence times measured by surface plasmon resonance (SPR), providing a fast comparative filter for lead optimization (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

### 4. Gaussian accelerated MD (GaMD) and LiGaMD
Miao and colleagues extended GaMD to Ligand Gaussian Accelerated MD (LiGaMD) (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf):
- Imposes harmonic boost potentials to the ligand–protein non-bonded potential and conformational bonded potential without requiring pre-selected collective variables.
- Enables spontaneous, repetitive ligand binding and unbinding cycles within hundreds of nanoseconds, allowing simultaneous calculation of $k_{\text{on}}$, $k_{\text{off}}$, and binding free energy through cumulant expansion reweighting (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

### 5. Weighted ensemble (WE) and Markov state models (MSMs)
- **Weighted ensemble (WESTPA)**: Divides reaction progress into geometric bins, running replicas in parallel and splitting/merging trajectories to maintain uniform occupancy. Captures rare unbinding transitions and probability flux without modifying the molecular Hamiltonian.
- **Markov state models (MSMs)**: Discretizes hundreds of microsecond-scale trajectories into conformational microstates, constructing transition probability matrices to compute mean first passage times (MFPT) and unbinding rates (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

## Benchmark target validations

The authors review validation benchmarks across model and therapeutic systems (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf):
- **T4 Lysozyme L99A**: Classic model for hydrophobic core binding, demonstrating that InMetaD accurately captures benzene egress pathways.
- **Kinases (CDK8/Cyclin C, p38 MAPK)**: Estimating prolonged residence times of Type-I vs Type-II kinase inhibitors associated with activation loop and DFG-out rearrangements.
- **GPCRs ($\beta_2\text{AR}, A_{2A}\text{AR}$)**: Characterizing unbinding pathways through extracellular vestibules and membrane lipid boundaries (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

## Related pages
- [[binding-kinetics]]
- [[enhanced-sampling]]
- [[molecular-dynamics-simulation]]
- [[binding-thermodynamics]]
- [[binding-free-energy]]
- [[protein-ligand-interactions]]

