# Binding Kinetics

**Summary**: Binding kinetics describes the rates at which a protein–ligand complex forms and dissociates; k_on and k_off determine K_d, and the dissociation rate — or its reciprocal, the residence time — is increasingly recognized as a predictor of in vivo efficacy.
**Sources**: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf; Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf
**Last updated**: 2026-11-01

---

## Basic kinetics

For P + L ⇌ PL (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

$$\frac{d[PL]}{dt} = k_{on}[P][L] - k_{off}[PL]$$

At equilibrium, k_on[P][L] = k_off[PL], so:

$$K_b = \frac{k_{on}}{k_{off}} = \frac{1}{K_d}$$

Units are M⁻¹·s⁻¹ for k_on and s⁻¹ for k_off. Fast association plus slow dissociation gives a high binding constant and high affinity. The residence time is τ = 1/k_off.

## Why kinetics matters

Binding free energy correlates with efficacy, but in vivo efficacy under nonequilibrium conditions is often better correlated with unbinding kinetics and residence time — the time a drug remains in contact with its target (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

Two drugs with the same K_d can have very different k_off: a compound with fast k_on and fast k_off equilibrates quickly, whereas a compound with slow k_on and slow k_off can remain bound long after free concentration falls. The latter can sustain target engagement between doses, but slow association may also delay onset (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### The translational residence time paradigm
The clinical impact of residence time ($t_R = 1/k_{\text{off}}$) is exemplified by human muscarinic $M_3$ receptor antagonists (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf):
- **Tiotropium** and **Ipratropium** share nearly identical picomolar binding affinities ($K_d \approx 10-100\text{ pM}$).
- However, tiotropium exhibits an exceptionally prolonged dissociation half-life ($t_R \approx 35\text{ hours}$), whereas ipratropium dissociates in approximately $0.3\text{ hours}$.
- This prolonged residence time permits tiotropium to be administered as a once-daily maintenance therapy for chronic obstructive pulmonary disease (COPD), whereas ipratropium requires multiple daily doses (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).
- In multi-kinase drug discovery (e.g., dasatinib, lapatinib), drugs engineered for slow off-rates against the therapeutic target and fast off-rates against off-target homologues maintain clinical efficacy while sharply reducing adverse toxicities (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf; Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

## Experimental measurement

### Surface plasmon resonance

SPR measures association and dissociation in real time. Protein is immobilized on a gold sensor surface; ligand flows over it; binding changes the refractive index in resonance units. Fitting the association and dissociation phases gives k_on and k_off, and K_b = k_on/k_off (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Limitations: mass transport limits accurate k_on measurement above ~10⁶ M⁻¹·s⁻¹, and immobilization perturbs the conformational and translational/rotational entropy of the protein, affecting the association rate (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Other methods

- **ITC** gives K_b and ΔH, hence ΔG and ΔS, but not rates directly.
- **Fluorescence polarization** gives equilibrium affinities via competition and can report kinetics in some formats.
- **Stopped-flow fluorescence** and **temperature-jump** methods resolve fast association and relaxation.
- **Single-molecule fluorescence** can observe individual binding events.

## Computational kinetics

### Unbiased MD

Unbiased MD is expensive but has important advantages: it requires no prior knowledge of binding poses, covers orthosteric and allosteric sites impartially, and provides detailed pathways and intermediate states (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

**Dasatinib and PP1 binding to Src kinase** (D. E. Shaw group, Anton):

- Dasatinib bound after 2.3 µs in one of four simulations (35 µs total) (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- PP1 bound after 15.1, 1.9, and 0.6 µs in seven simulations (115 µs total) (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- Dasatinib k_on = 4.3 s⁻¹·µM⁻¹, in striking agreement with the experimental 4.3 s⁻¹·µM⁻¹, although from a single event (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- Most successful binding events reproduced the crystallographic pose, supporting the funnel picture of binding (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- PP1 also docked at alternative pockets, showing that long unbiased runs can find allosteric sites (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- Water played a major role: crystallographic water positions were reproduced, and a water shell around PP1 created an entrance barrier (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

(source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md)

**β1/β2 adrenergic receptors:** 82 simulations, 230 µs total, 21 successful binding events. Alprenolol and dihydroalprenolol gave k_on = 31 s⁻¹·µM⁻¹, close to experiment. In one case FEP gave −13.4 ± 1.6 kcal/mol versus −12.2 kcal/mol experimentally. Of 12 bound poses, 6 matched the crystal structure; the others were less favorable prebinding poses. The simulations identified an **extracellular vestibule** as a major intermediate, and showed that the highest barrier was from bulk solvent to the vestibule, not from vestibule to site — caused by near-complete desolvation of the ligand and partial desolvation of the pocket (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Markov state models

MSMs convert many short trajectories into a transition matrix and yield rates, populations, and pathways (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

**Trypsin–benzamidine:** 495 simulations of 100 ns each; spontaneous binding in 187 cases (37%), with bound poses within 2 Å RMSD of the crystal structure. A 3D MSM with an 18 × 18 × 30 grid gave k_on and k_off in fair but not excellent agreement with experiment. Unbinding was never observed directly; k_off came from the MSM's equilibrium relations (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

**Apo-trypsin plasticity:** a 150 µs MSM identified six apo states, seven bound states, and four associated conformations. The six apo states were visited on tens-of-microseconds timescales and exposed different surface pockets, which were recognizable in crystal structures of wild-type and mutant protein (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Weighted ensemble / WExplore

WE runs many trajectories in parallel, pruning and spawning them across bins to maintain steady-state flux. It gives unbiased kinetics on timescales longer than the aggregated simulation (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

- **TPPU from soluble epoxide hydrolase:** experimental residence time 11 min (660 s); WExplore estimated 42 s with a standard error of ~10² s from only 6 µs of aggregated simulation — within 1–2 orders of magnitude. The analysis identified the interactions that most affect unbinding and the role of hydrogen bonding and solvation (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Benzamidine from trypsin:** predicted exit time 180 µs versus measured 1700 µs, and revealed three exit channels, one previously unknown, two involving large loop motions (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **FK506 ligands:** nanosecond residence times, where WExplore matched plain MD at much lower cost and mapped ligand exit-point distributions (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Metadynamics and Infrequent Metadynamics (InMetaD)

Standard metadynamics accelerates unbinding along chosen collective variables (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Infrequent metadynamics (InMetaD)**: Tiwary and Parrinello established that by depositing repulsive Gaussian bias hills at a frequency significantly slower than the transition-state crossing time ($\tau_{\text{dep}} \gg \tau_{\text{trans}}$), bias accumulates only in the bound basin without corrupting the transition barrier (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).
- **Poisson statistics**: Unbinding times across independent replicas follow a Poisson distribution, validated by the Kolmogorov–Smirnov test ($p > 0.05$). The physical unbinding rate is recovered by multiplying the biased transition time by the time-averaged acceleration factor $\alpha = \langle \exp(V_{\text{bias}}/k_B T) \rangle$ (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).
- **OPES-flooding**: An alternative that rapidly and stably floods the bound basin to accelerate transition-state escape with minimal tuning (source: Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

One HSP90 study used 150–750 ns runs and estimated a residence time of **21 s**, in excellent agreement with experiment, for a system where 35 µs of unbiased MD never observed unbinding (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Ligand Gaussian accelerated MD (LiGaMD)

Developed by Miao and coworkers, Ligand Gaussian Accelerated Molecular Dynamics applies harmonic boost potentials to ligand–protein non-bonded interaction potentials and torsional dihedrals without requiring pre-defined collective variables (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf; Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf):
- Facilitates repetitive, spontaneous ligand association and dissociation events within hundreds of nanoseconds of simulation time (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).
- By reweighting trajectories via 2nd-order cumulant expansion, LiGaMD computes both binding free energies ($\Delta G_{\text{bind}}$) and kinetic rate constants ($k_{\text{on}}$, $k_{\text{off}}$, $\tau = 1/k_{\text{off}}$) simultaneously within a single simulation framework (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

### Supervised molecular dynamics (SuMD)

Developed by Sabbadin and Moro, SuMD monitors ligand–pocket distance across consecutive short classical MD windows (e.g., 600 ps) (source: Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf):
- If the ligand approaches the binding site, the simulation continues; if not, velocities are randomly reassigned from the previous state.
- Reconstructs complete association pathways, transient intermediate vestibules, and desolvation events within nanoseconds without adding energetic biases to the Hamiltonian (source: Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf).

### Scaled MD

Scaled MD reduces the potential by a factor λ, flattening barriers (see [[enhanced-sampling]]). With λ = 0.4 — equivalent to a 2.5-fold temperature increase — all ligands of HSP90, Grp78, and the A2A GPCR unbinding within nanoseconds to tens of nanoseconds, with backbone restraints outside the site preventing unfolding. Repeats give an average unbinding time used to rank compounds (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

For seven noncongeneric glucokinase activators, λ = 0.5 was adequate in all cases, while λ = 0.4 was sufficient for congeneric ranking (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### τ-RAMD

τ-RAMD applies a random force to the ligand and reassigns its direction when progress stalls. It accelerates unbinding by orders of magnitude; 40–200 simulations per compound give a mean unbinding time τ. The average deviation from experimental residence times was **2.3τ**, reduced to **2τ** for congeneric compounds (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

### Other approaches

- **ABMD with an electrostatic collective variable** has been used for residence-time ranking (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Diffusion models** and **Brownian dynamics** estimate diffusion-limited association rates and encounter pathways (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Machine learning** on MD features is increasingly used to classify binding/unbinding pathways (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## What determines k_off

From the simulations above, recurrent factors are (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):

- **Desolvation** — a ligand must shed water to enter or leave a buried site; the bulk-to-vestibule barrier in β-AR was dominated by ligand desolvation.
- **Water networks** — bridging waters can either stabilize the bound state or create an exit barrier.
- **Loop and side-chain gating** — exit channels often open through large-amplitude loop motions.
- **Multiple metastable states** — intermediate stations can slow dissociation without appearing in a single structure.
- **Electrostatics** — long-range attraction can speed association; charge burial can slow dissociation.

## Practical guidance

- Measure both k_on and k_off when possible; a single K_d hides kinetic differences.
- For slow-dissociating compounds, use enhanced sampling (WE, metadynamics, scaled MD, τ-RAMD) rather than plain MD.
- Validate computational rankings on a congeneric series, where relative errors are smaller.
- Report the method, collective variables, number of replicas, and force field; kinetic predictions are method-dependent.
- Treat absolute rates from accelerated methods as order-of-magnitude estimates unless validated.

## Related pages
- [[binding-thermodynamics]]
- [[binding-free-energy]]
- [[experimental-binding-assays]]
- [[benchmark-drug-targets]]
- [[enhanced-sampling]]
- [[molecular-dynamics-simulation]]
- [[water-and-solvation]]
- [[decherchi-2020-binding-thermodynamics-kinetics]]
- [[du-2016-protein-ligand-interactions]]
- [[ahmad-2022-enhanced-sampling-binding-kinetics]]
- [[adediwura-2024-free-energy-kinetics-drug-discovery]]
- [[metadynamics]]
- [[markov-state-models]]
