# Decherchi & Cavalli (2020) — Thermodynamics and Kinetics of Drug–Target Binding

**Summary**: A comprehensive review of the statistical mechanics, enhanced-sampling methods, and free-energy and kinetics calculations used to predict drug–target binding, with an emphasis on which methods are mature enough for prospective drug discovery.
**Sources**: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md
**Last updated**: 2026-11-01

---

## Citation and scope

Sergio Decherchi and Andrea Cavalli, "Thermodynamics and Kinetics of Drug–Target Binding by Molecular Simulation," *Chemical Reviews* 120(23):12788–12833 (2020), doi:10.1021/acs.chemrev.0c00534. PMC8011912 (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

The review covers the theory and practice of simulating drug–target binding at atomic resolution. Its scope is deliberately limited: it treats the drug–target complex in solution or in a membrane and neglects organism-level pharmacology. It is the wiki's most detailed source for [[enhanced-sampling]], [[binding-free-energy]], and [[binding-kinetics]].

## The problem

A complete description must span femtoseconds (molecular vibrations) to milliseconds and beyond (binding and unbinding). With a 1 fs time step, a millisecond event requires 10¹² integration steps and a second requires 10¹⁵ — far beyond current hardware for routine work. The difficulty is not only raw speed but near-nonergodicity: the free-energy landscape is divided into basins separated by high barriers, and the system may never cross them in an affordable simulation (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

The review separates two goals that are often conflated (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):

1. **Equilibrium thermodynamics** — sample the Boltzmann distribution and compute free-energy differences.
2. **Kinetics** — quantify rates under stationary or nonequilibrium conditions, which is harder.

## Structure of the review

1. Introduction and scope.
2. Basics of molecular simulation: models, force fields, MD.
3. Enhanced sampling: steered MD, ABMD, parallel tempering, scaled MD, τ-RAMD, metadynamics, MSMs, weighted ensemble.
4. Free-energy methods: umbrella sampling, ABF, relative FEP/TI, double annihilation/decoupling, MM-PBSA/GBSA.
5. Applications: absolute binding free energy, relative binding free energy, binding and unbinding kinetics.
6. Recent machine-learning trends.
7. Practical guidelines.
8. Conclusions and perspectives.

## Enhanced sampling

### Steered MD and ABMD

Steered MD adds a harmonic potential whose center moves along the reaction coordinate at constant velocity, pulling the ligand out of the site. Park and Schulten showed how to recover the potential of mean force from several nonequilibrium replicas using Jarzynski's equality (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Adiabatic bias MD (ABMD) moves a harmonic restraint only when the collective variable advances in the desired direction — a pawl-and-ratchet mechanism. Its speed is set partly by the system's natural evolution, so it produces gentler, more physically plausible trajectories than steered MD (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Parallel tempering

Replicas at different temperatures exchange configurations with a Metropolis acceptance probability. Equal acceptance requires a geometric temperature progression. The method satisfies detailed balance and can also exchange other observables, such as metadynamics bias potentials (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Scaled MD

A positive constant μ ∈ (0, 1] defines U_μ(r) = μU(r), uniformly flattening all barriers. Constant scaling is equivalent to a 1/μ temperature scaling. Backbone restraints outside the binding site prevent unfolding. The unbinding time is defined as the time at which the ligand becomes fully surrounded by water; at least 20 repeats define an average. The method needs few free parameters and no reaction coordinate, but requires prior knowledge of the binding-site residues and can distort the pathway under heavy scaling (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

A published k_off ranking used **λ = 0.4** — equivalent to a 2.5-fold temperature increase — on HSP90, Grp78, and the A2A GPCR with three ligand series each. For seven noncongeneric glucokinase activators, λ = 0.5 was adequate; λ = 0.4 sufficed for congeneric ranking (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### τ-RAMD

Random acceleration MD periodically applies a random force to the ligand and reassigns its direction when progress stalls. τ-RAMD runs on average **40–200 simulations per compound** and uses the mean unbinding time τ to correlate with experimental residence times. Average deviation was **2.3τ**, reduced to **2τ** for congeneric compounds (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Metadynamics

Gaussians are deposited along collective variables; the negative of the converged bias approximates the free-energy surface. The method explores CV space automatically and yields free energy at convergence, but it is hard to know when to stop and the last Gaussian introduces an offset. **Well-tempered metadynamics** reduces the Gaussian height over time and has been proved to converge (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Markov state models and weighted ensemble

An MSM partitions configuration space into states and estimates a transition matrix over a lag time; the Chapman–Kolmogorov test validates the lag. MSMs give populations, free energies, and rates, but work best for diffusive problems and poorly for high-barrier activated events (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Weighted ensemble (WE) partitions the path into bins, runs many trajectories, and prunes/spawns walkers to maintain steady-state flux. It provides unbiased kinetics on timescales longer than the aggregated simulation and does not require a Markovian assumption. **WExplore** defines bins dynamically as hierarchical Voronoi polyhedra (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Free-energy methods

### Umbrella sampling and ABF

Umbrella sampling adds a harmonic restraint at successive positions along the reaction coordinate; windows are recombined with WHAM or MBAR. To compare with experiment, the PMF must be converted to a standard binding free energy, correcting from the standard-state volume **V₀ = 1661 Å³** (1 M) to the unbound volume sampled and accounting for restraints orthogonal to the coordinate (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Adaptive biasing force (ABF) estimates the mean force on the coordinate on the fly and applies an opposing bias; the free energy follows by integration (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Relative FEP and TI

A thermodynamic cycle gives:

$$\Delta\Delta G_{bind} = \Delta G_{site}(a\to b) - \Delta G_{solv}(a\to b)$$

Single topology mutates atoms in place; dual topology keeps both end states and interpolates U(λ) = (1−λ)U_a + λU_b. Soft-core potentials with parameters α_LJ and α_C avoid endpoint clashes (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

FEP/REST restricts the "hot" region to the ligand and nearby residues, greatly reducing the number of replicas compared with parallel tempering. A retrospective assessment on **200 ligands and 10 targets** performed **330 perturbations**, many changing up to 10 non-hydrogen atoms, with an absolute error below **1.5 kcal/mol in 81.2%** of cases (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

TI integrates the ensemble average of ∂U/∂λ over λ windows. It is exact in the limit of infinite sampling and a good force field, and converges with simulation time — unlike MM-GBSA, whose convergence is harder to assess (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Absolute binding free energy

**Double annihilation** annihilates the ligand in water and in the complex. It ignores the standard-state correction, and the decoupled ligand can wander through the whole box, making convergence difficult.

**Double decoupling** adds restraints that keep the ligand in the site during decoupling, yielding a rigorous standard binding free energy. The restraint contribution must be computed; an analytical estimate due to Karplus and coworkers makes the scheme practical (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Absolute alchemical calculations reach about **1–2 kcal/mol** accuracy versus ITC, at a cost of 10³–10⁴ core hours per cycle (2019) (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### MM-PBSA and MM-GBSA

These estimate:

$$\Delta G_{bind} = G_{complex} - G_{protein} - G_{ligand}, \qquad G \approx E_{MM} + G_{pol} + G_{np} - TS$$

The single-trajectory variant reuses one trajectory for all three species, roughly halving cost and exploiting error cancellation; it is valid only if binding causes no major conformational change. The entropy term is the hardest part — normal-mode or quasiharmonic analysis scales poorly, and ΔS is often set to zero (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Accuracy: MM-PBSA could not consistently estimate absolute binding free energy within an RMSD of about **5 kcal/mol** on bromodomains; with entropy estimation it recovered about **90% of the ranking accuracy** of alchemical methods at **5% of the cost**. Results depend on force field, charges, PB versus GB, single versus multiple trajectories, and explicit waters. PB is theoretically superior to GB, but factor Xa studies found the opposite. Single-snapshot predictions sometimes outperform trajectory averages, and ~200 ps runs sometimes outperform ~500 ps runs — signs of incomplete convergence (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Applications

### Absolute binding free energy

A bromodomain study estimated absolute binding free energies for 11 inhibitors. With experimental geometries the error was about **1 kcal/mol**; starting from docking poses it increased only slightly and remained close to 1 kcal/mol, far better than docking alone. Each cycle required microsecond simulations and 10³–10⁴ core hours, excluding extensive screening but making alchemical methods attractive for rescoring (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

A 36-complex study across 22 bromodomains with three ligands achieved a standard deviation below 1 kcal/mol and a maximum deviation below 2 kcal/mol for the RVX ligands; bromosporine was less accurate, with a standard error of about 2 kcal/mol (a factor of 30 in concentration), and improved when AM1-BCC charges were replaced by RESP and torsions were refined (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Kinetics from unbiased MD

**Dasatinib and PP1 binding to Src kinase** (D. E. Shaw group, Anton):

- Dasatinib bound after 2.3 µs in one of four simulations (35 µs total).
- PP1 bound after 15.1, 1.9, and 0.6 µs in seven simulations (115 µs total).
- Dasatinib k_on = **4.3 s⁻¹·µM⁻¹**, in striking agreement with the experimental 4.3 s⁻¹·µM⁻¹, although from a single event.
- Most successful events reproduced the crystallographic pose.
- PP1 also bound at alternative pockets, showing that long unbiased runs can find allosteric sites.
- Water played a major role: crystallographic water positions were reproduced, and a water shell around PP1 created an entrance barrier.

(source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md)

**β1/β2 adrenergic receptors:** **82 simulations, 230 µs total, 21 successful binding events**. Alprenolol and dihydroalprenolol gave k_on = 31 s⁻¹·µM⁻¹, close to experiment; in one case FEP gave −13.4 ± 1.6 kcal/mol versus −12.2 kcal/mol experimentally. Of 12 bound poses, 6 matched the crystal structure. The simulations identified an **extracellular vestibule** as a major intermediate and showed that the highest barrier was from bulk solvent to the vestibule, caused by near-complete ligand desolvation and partial pocket desolvation (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Kinetics from MSMs and WE

**Trypsin–benzamidine:** 495 simulations of 100 ns each; spontaneous binding in 187 cases (**37%**), with bound poses within 2 Å RMSD. A 3D MSM with an 18 × 18 × 30 grid gave k_on and k_off in fair agreement with experiment; unbinding was never observed directly, and k_off came from equilibrium relations (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

**Apo-trypsin plasticity:** a 150 µs MSM identified six apo states, seven bound states, and four associated conformations; the six apo states exposed different surface pockets recognizable in crystal structures of wild-type and mutant protein (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

**TPPU from soluble epoxide hydrolase:** experimental residence time 11 min (**660 s**); WExplore estimated **42 s** with a standard error of ~10² s from only 6 µs of aggregated simulation — within 1–2 orders of magnitude (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

**Benzamidine from trypsin:** predicted exit time 180 µs versus measured 1700 µs, revealing three exit channels, one previously unknown (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Metadynamics kinetics

An HSP90 study used 150–750 ns runs and estimated a residence time of **21 s**, in excellent agreement with experiment, for a system where 35 µs of unbiased MD never observed unbinding (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Machine-learning trends

The review briefly covers deep learning for affinity prediction and ML approaches for learning collective variables and Boltzmann distributions. It treats these as promising but not yet a replacement for physics-based free-energy calculations (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[machine-learning-scoring-functions]].

## Practical guidelines

The review's most valuable section for practitioners (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):

1. **Start from solid structural data.** High-resolution crystal or co-crystal structures, preferably better than 2.5 Å, were considered essential by the authors, though subsequent benchmarks (Muegge & Hu 2023, van Pinxteren & Jespers 2025) demonstrate that homology models and AlphaFold holo complexes relaxed via equilibrium MD achieve comparable quantitative accuracy.
2. **Prefer methods with few user choices.** If a method requires many human decisions that significantly change the outcome, it is not mature enough for prospective quantitative use.
3. **Beginners should start with methods that need no collective variable** — FEP and scaled MD — before moving to umbrella sampling or metadynamics.
4. **Match method to budget and prior knowledge.**
   - Large budget and good prior knowledge: umbrella sampling, metadynamics, MSMs.
   - Time-constrained: FEP, scaled MD.
   - Very limited resources or virtual screening: scoring functions or MM-PBSA/GBSA.
5. **Exploration:** metadynamics is the best choice when a reaction coordinate can be guessed; parallel tempering is best for pure exploration without a CV; WE is generally preferable to MSM and can replace it in almost all cases.
6. **Absolute free energy:** double decoupling is the rigorous choice because it needs no CV, but restraints and cost require care. Metadynamics and umbrella sampling can also give absolute values, but converting a PMF into a binding free energy requires partitioning bound and unbound regions, which is ambiguous and user-dependent.
7. **Relative free energy:** FEP and TI are mature enough for systematic industrial use. Their limitations are large perturbations (scaffold hopping) and net charge changes.
8. **Kinetics:** MSM and WE estimate rates directly; scaled MD and τ-RAMD are much faster for k_off ranking. In congeneric series, k_on is often nearly constant, so k_off ranking correlates with K_D and free energy.
9. **Always run replicates** and report convergence, because a single trajectory can differ from another by as much as the full range of interest.

## Limitations and open problems

- **Prospective validation is rare.** Most methods are applied retrospectively without blind testing; standardized benchmarks, blind competitions, common force fields, shared code, and data repositories are needed (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Force-field divergence.** Results can differ between force fields, and there is no common open gold standard (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Sampling and computational power.** One simulation is not enough to assess error; many replicates are needed, and current hardware is insufficient for robust routine work (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Standardization.** Even residue naming is not standardized (HIS, HIP, HI+ for double-protonated histidine), a symptom of a broader problem (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **MM-PBSA/GBSA convergence.** The method is not systematically convergent, and its entropy treatment is crude (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **PMF-to-ΔG conversion.** Partitioning bound and unbound states is ambiguous and can dominate the result (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Key numbers

| Quantity | Value |
|---|---|
| 1 ms at 1 fs | 10¹² steps |
| Standard-state volume | 1661 Å³ (1 M) |
| FEP/REST assessment | 200 ligands, 10 targets, 330 perturbations |
| FEP/REST accuracy | <1.5 kcal/mol in 81.2% |
| Absolute alchemical accuracy | ~1–2 kcal/mol vs ITC |
| MM-PBSA absolute error | ~5 kcal/mol on bromodomains |
| MM-PBSA ranking | ~90% of alchemical accuracy at 5% cost |
| Dasatinib k_on | 4.3 s⁻¹·µM⁻¹ (experiment 4.3) |
| β-AR simulations | 82 runs, 230 µs, 21 binding events |
| TPPU residence time | 660 s experimental; 42 s WExplore |
| τ-RAMD runs | 40–200 per compound |
| τ-RAMD deviation | 2.3τ (2τ congeneric) |
| Scaled MD λ | 0.4 |

(source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md)

## Strengths

- Unusually thorough coverage of both thermodynamics and kinetics, with theory and application.
- Gives quantitative accuracy and cost for each method, which is rare in reviews.
- Explicit, opinionated practical guidelines — including which methods are not yet mature.
- Honest about retrospective bias, force-field divergence, and standardization failures.

## Limitations

- Dense and technical; readers without statistical mechanics will struggle with the theory sections.
- The applications are selective and organized by method rather than by target class.
- Machine learning is covered only briefly.
- Published in 2020; hardware, ML potentials, and enhanced-sampling software have advanced since.

## Why it matters for this wiki

This is the main source for [[enhanced-sampling]], [[binding-free-energy]], and [[binding-kinetics]], and a supporting source for [[molecular-dynamics-simulation]], [[force-fields]], [[water-and-solvation]], and [[cryptic-and-allosteric-pockets]]. It supplies the wiki's accuracy hierarchy — docking < MM-PBSA/GBSA < FEP/TI — and the concrete benchmarks (81.2% within 1.5 kcal/mol, 4.3 s⁻¹·µM⁻¹, 42 s versus 660 s) that make the methods comparable.

## Related pages
- [[enhanced-sampling]]
- [[binding-free-energy]]
- [[binding-kinetics]]
- [[molecular-dynamics-simulation]]
- [[force-fields]]
- [[water-and-solvation]]
- [[binding-thermodynamics]]
- [[du-2016-protein-ligand-interactions]]
- [[hollingsworth-2018-md-for-all]]
- [[durrant-2011-md-drug-discovery]]
