# Adediwura et al. (2024) — Binding Free Energy and Kinetics in Modern Drug Discovery

**Summary**: An expert review from the Miao laboratory unifying alchemical free energy, MM-PBSA/GBSA, and Gaussian Accelerated MD kinetics calculations within modern pharmacology.
**Sources**: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Victor A. Adediwura, Kushal Koirala, Hung N. Do, Jinan Wang, and Yinglong Miao, "Understanding the impact of binding free energy and kinetics calculations in modern drug discovery," *Expert Opinion on Drug Discovery* 19(5):543–556 (2024), doi:10.1080/17460441.2024.2334709; author manuscript PMC11090448 / NIHMS1990962 (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

The review synthesizes computational techniques for calculating drug binding thermodynamics ($\Delta G_{\text{bind}}$) and drug binding kinetics ($k_{\text{on}}$, $k_{\text{off}}$, and residence time $t_R = 1/k_{\text{off}}$). Originating from the developers of Gaussian Accelerated Molecular Dynamics (GaMD and LiGaMD), the authors demonstrate how combining thermodynamic affinity with kinetic rate prediction guides hit-to-lead progression, target selectivity, and translational in vivo efficacy (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

## The computational binding hierarchy

The authors organize computational binding methodologies into an overarching four-tier hierarchy based on accuracy, physical rigor, and computational cost (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf):

### 1. Alchemical free energy calculations (FEP and TI)
Alchemical approaches—Free Energy Perturbation (FEP) and Thermodynamic Integration (TI)—represent the gold standard for congeneric lead optimization:
- Employ unphysical $\lambda$-paths and soft-core potentials to bypass endpoint singularities at $\lambda = 0$ or $1$.
- Evaluate relative binding free energies ($\Delta\Delta G_{\text{bind}}$) with prospective errors $\le 1.0\text{ kcal/mol}$.
- Highly demanding computationally, restricting routine deployment to dozens or hundreds of congeneric molecules (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

### 2. Endpoint free energy methods (MM-PBSA and MM-GBSA)
Endpoint methods offer a strategic compromise between speed and physical accuracy:
- Partition binding free energy into gas-phase molecular mechanics ($E_{\text{MM}}$), continuum polar solvation (Poisson–Boltzmann or Generalized Born), non-polar surface terms (SASA), and conformational entropy.
- **Single-trajectory approach**: Extracting receptor and ligand configurations from a single MD complex trajectory eliminates huge fluctuations in internal bonded energies ($\Delta E_{\text{bonded}} = 0$), yielding higher precision.
- **Entropy treatment**: Normal Mode Analysis (NMA) is computationally burdensome and introduces noise. While Interaction Entropy (IE) enhances correlation in certain targets (such as Bcl-2 family members), other benchmarks (such as wild-type and mutant EGFR) indicate that IE can introduce noise, confirming that entropy parameter tuning remains target-dependent (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

### 3. Physical pathway methods (US, SMD, ABF)
Umbrella Sampling (US), Steered MD (SMD), and Adaptive Biasing Force (ABF) track the literal spatial entry and exit of the drug along a physical reaction coordinate, constructing Potential of Mean Force (PMF) profiles that reveal free-energy barriers and intermediate binding vestibules (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

### 4. Advanced kinetic simulations (GaMD, LiGaMD, InMetaD, MSMs)
- **Gaussian Accelerated MD (GaMD)**: Adds harmonic boost potentials that follow near-Gaussian distributions, enabling unbiasing via cumulant expansion.
- **Ligand GaMD (LiGaMD)**: Applies targeted boost potentials to both the ligand–protein non-bonded interaction potential and ligand bonded/torsional potential. LiGaMD captures repetitive, spontaneous ligand binding and dissociation events on microsecond timescales, permitting simultaneous calculation of binding free energies ($\Delta G_{\text{bind}}$), association rate constants ($k_{\text{on}}$), and dissociation rate constants ($k_{\text{off}}$) without requiring predefined reaction coordinates (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).
- **Other kinetic engines**: Infrequent Metadynamics (InMetaD), $\tau$-Random Acceleration MD ($\tau$-RAMD), Weighted Ensemble (WESTPA), and Markov State Models (MSMs) (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

## Translational and pharmacological impact of residence time

The authors detail why modern drug discovery increasingly prioritizes drug residence time ($t_R = 1/k_{\text{off}}$) alongside equilibrium affinity ($K_d$) (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf):
- **Pharmacodynamics under clearance**: In vivo drug concentrations decline continuously following dosing. A drug with a slow off-rate remains bound to its target long after free drug has been eliminated from circulation, sustaining target engagement.
- **Clinical paradigm — Muscarinic antagonists**:
  - *Tiotropium* and *Ipratropium* possess nearly identical equilibrium binding affinities for the human muscarinic $M_3$ receptor ($K_d \approx 10-100\text{ pM}$).
  - However, tiotropium has an exceptionally prolonged dissociation half-life ($t_R \approx 35\text{ hours}$), whereas ipratropium dissociates rapidly ($t_R \approx 0.3\text{ hours}$).
  - This kinetic distinction allows tiotropium to function as a once-daily maintenance therapy for chronic obstructive pulmonary disease (COPD), while ipratropium requires frequent daily dosing (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).
- **Target selectivity and toxicity reduction**: In multi-target families such as protein kinases, designing drugs (e.g., lapatinib, dasatinib) with slow dissociation kinetics against the disease target and fast dissociation kinetics against off-target homologues minimizes systemic adverse effects (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

## Related pages
- [[binding-kinetics]]
- [[binding-free-energy]]
- [[enhanced-sampling]]
- [[molecular-dynamics-simulation]]
- [[binding-thermodynamics]]
- [[structure-based-drug-design]]

