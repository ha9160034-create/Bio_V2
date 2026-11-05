# Binding Free Energy

**Summary**: Binding free energy (ΔG_bind) quantifies the strength of a protein–ligand complex and is related to the binding constant by ΔG° = −RT ln K_b; it can be computed by alchemical free-energy perturbation or thermodynamic integration, by potential-of-mean-force methods, or approximately by MM-PBSA/MM-GBSA.
**Sources**: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulations and drug discovery.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Molecular dynamics simulation for all.md; Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf; Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf; Integrating Machine Learning into Free Energy Perturbation Workflows.pdf; The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf; Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf; The Open Force Field Initiative Open Software and Open Science.pdf; Training a force field for proteins and small molecules from.pdf; Best Practices for Foundations in Molecular Simulations [Article.pdf
**Last updated**: 2026-11-02

---

## Definition and relationship to experiment

The standard binding free energy is the free-energy change at 1 atm, 298 K, and 1 M effective reactant concentrations. It is related to the binding constant by (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

$$\Delta G^\circ = -RT\ln K_b$$

with R = 1.987 cal·K⁻¹·mol⁻¹. A larger K_b gives a more negative ΔG°. Since K_b = k_on/k_off = 1/K_d, kinetics and thermodynamics are linked.

Free energy is a state function: ΔG depends only on the initial and final states, not on the path. This is what licenses alchemical transformations that follow nonphysical paths (source: Molecular dynamics simulations and drug discovery.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Why compute it

Docking scores are fast but approximate. Free-energy calculations aim to predict affinity accurately enough to guide synthesis. The practical hierarchy, from cheap to expensive (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):

1. docking scoring functions (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md);
2. MM-PBSA / MM-GBSA (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md);
3. alchemical FEP/TI and potential-of-mean-force methods (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Relative binding free energies are enough for lead optimization within a congeneric series. Absolute binding free energies are needed to compare chemically distinct scaffolds, to assess selectivity across proteins, and to estimate intrinsic potency (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Alchemical relative free energy: FEP

For two ligands a and b, a thermodynamic cycle relates the relative binding free energy to the alchemical transformations in the binding site and in solution:

$$\Delta\Delta G_{bind} = \Delta G_{site}(a\to b) - \Delta G_{solv}(a\to b)$$

The horizontal (binding) legs are hard; the vertical (mutation) legs are feasible perturbations, which is the basis of FEP for ligand series (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Two topology schemes are used:

- **Single topology** — each atom mutates in place; force-field parameters are interpolated between endpoints. It minimizes the number of transformations, but bonds are elongated and shortened, which is not a weak perturbation.
- **Dual topology** — both end states coexist in scaled form, and the potential is an interpolation U(λ) = (1−λ)U_a + λU_b. Only endpoint parameters are needed, but emerging atoms can clash with vanishing atoms near the endpoints.

The clash problem is handled by a **soft-core potential** that modifies the Lennard-Jones term with a positive constant α_LJ and an analogous modification of the Coulomb term with α_C (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### FEP/REST

FEP with replica exchange with solute tempering (REST) restricts the "hot" region to the ligand and nearby residues, excluding bulk water. This greatly reduces the number of replicas compared with parallel tempering and improves sampling of dihedral barriers. Each λ window is assigned a solute temperature, rising to T_h at the middle window and returning to T₀ at the end (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

A retrospective FEP/REST assessment on 200 ligands and 10 targets performed **330 perturbations**, many changing up to 10 non-hydrogen atoms, with an absolute error below 1.5 kcal/mol in **81.2%** of cases (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Alchemical absolute free energy

### Double annihilation

The ligand is annihilated in water (ΔG₁) and in the complex (ΔG₂), giving ΔG_bind = ΔG₂ − ΔG₁. The method has two problems: it ignores the standard-state correction, so results are not strictly comparable with experiment, and after decoupling the ligand can wander through the whole box, making convergence very difficult (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Double decoupling

Gilson's rigorous version adds restraints that keep the ligand in the binding site during decoupling, avoiding the sampling problem and yielding a correct standard binding free energy. The restraint free-energy contribution must be computed; for protein–ligand binding, an analytical estimate due to Karplus and coworkers makes the scheme practical (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Thermodynamic integration

TI computes the free-energy difference by integrating the ensemble average of the derivative of the Hamiltonian with respect to the coupling parameter:

$$\Delta G = \int_0^1 \left\langle \frac{\partial U(\lambda)}{\partial \lambda}\right\rangle_\lambda d\lambda$$

The integrand is evaluated at a series of λ windows. TI is exact in the limit of infinite sampling and a good force field; with increasing simulation time it converges to that limit, unlike MM-GBSA, whose convergence properties are harder to assess (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Benchmark Design and Force Field Advances in FEP
In establishing benchmark standards for binding free energy calculations, Hahn et al. (LiveCoMS, 2022) formulate essential criteria for constructing, preparing, and evaluating alchemical datasets (PLBenchmarks) (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf):
- **Direct Biophysical Assays**: Prioritize ITC and SPR measurements over functional IC50 data; IC50 conversions via Cheng–Prusoff are valid only when enzyme concentration is negligible compared to Ki, and variations in substrate conditions distort affinity comparisons (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Perturbation Graph Topology**: Alchemical networks must incorporate closed, redundant cycles; cycle closure hysteresis ($\sum_{\text{cycle}} \Delta\Delta G = 0$) provides an automated, intrinsic convergence test across perturbations (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Dynamic Range and Experimental Noise Floor**: Benchmarks require an experimental affinity span of $> 2\text{--}3$ orders of magnitude ($\Delta\Delta G \ge 2.7\text{--}4.1\text{ kcal/mol}$), noting that experimental measurement noise sets an unavoidable error floor of $\sim 0.3\text{--}0.5\text{ kcal/mol}$ (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Modern Open Force Fields**: OpenFF Sage 2.2.0 (employing direct chemical perception via SMIRNOFF) and Garnet (employing continuous graph neural networks with double exponential nonbonded potentials) achieve RBFE accuracy competitive with proprietary engines like OPLS4 (RMSE $\approx 1.1\text{--}1.4\text{ kcal/mol}$) across diverse pharmaceutical targets (source: The Open Force Field Initiative Open Software and Open Science.pdf; Training a force field for proteins and small molecules from.pdf). See [[open-force-field-and-ligand-parameterization]] and [[simulation-and-benchmark-best-practices]].

## Potential-of-mean-force methods

Umbrella sampling, metadynamics, steered MD, and related methods compute a free-energy profile along a physical coordinate. They have the advantage of providing kinetic information in addition to thermodynamics (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

The PMF must be converted to a standard binding free energy:

- integrate the PMF over the bound and unbound regions to obtain the ratio of partition functions (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md);
- correct from the standard-state volume V₀ = 1661 Å³ (1 M) to the unbound volume actually sampled (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md);
- account for restraints applied orthogonal to the reaction coordinate (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Without these corrections, PMF values are not comparable with experimental ΔG° (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[enhanced-sampling]].

## MM-PBSA and MM-GBSA

These methods estimate the binding free energy from a thermodynamic cycle:

$$\Delta G_{bind} = G_{complex} - G_{protein} - G_{ligand}$$

Each term is approximated as:

$$G \approx E_{MM} + G_{pol} + G_{np} - TS$$

where E_MM contains bonded, van der Waals, and electrostatic molecular-mechanics terms; G_pol is the polar solvation free energy from the Poisson–Boltzmann equation (MM-PBSA) or Generalized Born (MM-GBSA); G_np is a nonpolar term, often proportional to solvent-accessible surface area; and −TS is the entropic contribution (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Practical features (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):

- They are post-processing steps on an MD or MC trajectory with explicit solvent.
- The **single-trajectory** variant extracts complex, protein, and ligand configurations from the same trajectory, roughly halving the cost and exploiting error cancellation. It is valid only if binding causes no major conformational change; otherwise a multi-trajectory approach is advised.
- The **entropy term** is the hardest part. Normal-mode or quasiharmonic analysis requires building and diagonalizing a Hessian, which scales poorly with atom count. Often ΔS is set to zero, which is a drastic assumption.
- Water molecules in the pocket contribute an entropy change that is difficult to model and is a major source of error.

Accuracy and caveats:

- MM-PBSA is generally preferable to docking scoring functions, but it is questionable for choosing a substituent before synthesis (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- On bromodomains, MM-PBSA could not consistently estimate absolute binding free energy within an RMSD of about 5 kcal/mol in the best case; with entropy estimation it recovered about 90% of the ranking accuracy of alchemical methods at 5% of the cost (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- Results depend on force field, charges, PB versus GB, single versus multiple trajectories, and whether explicit waters are included. PB usually outperforms GB theoretically, but factor Xa studies found the opposite, illustrating the lack of systematic trends (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). Systematic benchmarks by Genheden & Ryde confirm that MM-GBSA often yields ranking correlations equal or slightly superior to MM-PBSA because analytical Generalized Born continuum electrostatics avoid numerical grid-discretization errors and benefit from error cancellation across trajectory snapshots (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf). See [[mm-pbsa-and-mm-gbsa]].
- Single-snapshot predictions sometimes outperform trajectory averages, and short ~200 ps runs sometimes outperform ~500 ps runs — a sign that the method is not systematically convergent (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Accuracy and cost

| Method | Typical accuracy | Cost | Best use |
|---|---|---|---|
| Docking score | qualitative; poor affinity correlation | seconds–minutes per ligand | ultra-large screening |
| MM-PBSA/GBSA | ~5 kcal/mol absolute; useful ranking | minutes–hours per ligand | rescoring, large series |
| Relative FEP/TI | ≤1.5 kcal/mol in 81.2% of cases | hours–days per perturbation | congeneric lead optimization |
| Absolute alchemical | ~1–2 kcal/mol versus ITC | 10³–10⁴ core hours per cycle (2019) | selectivity, distinct scaffolds |
| PMF/umbrella/metadynamics | ~1–2 kcal/mol with good sampling | high | mechanism plus thermodynamics |

Sources: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md.

A bromodomain absolute free-energy study illustrates the state of the art: with experimental geometries the error was about 1 kcal/mol; starting from docking poses the error increased only slightly and remained close to 1 kcal/mol, far better than docking alone. Each alchemical cycle required microsecond simulations and 10³–10⁴ core hours, which excludes extensive screening but makes alchemical methods attractive for rescoring (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

A 36-complex study across 22 bromodomains with three ligands achieved a standard deviation below 1 kcal/mol and a maximum deviation below 2 kcal/mol for the RVX ligands; bromosporine was less accurate, with a standard error of about 2 kcal/mol (a factor of 30 in concentration), and improved when AM1-BCC charges were replaced by RESP and torsions were refined (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Modern alchemical estimators and error control

Modern alchemical free energy calculations rely on rigorous minimum-variance estimators (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf):
- **Multistate Bennett Acceptance Ratio (MBAR)**: Outperforms simple FEP exponential averaging and TI numerical integration by incorporating configurations from all sampled $\lambda$-windows simultaneously. MBAR eliminates histogram binning and provides analytical asymptotic variances and phase space overlap matrices (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).
- **Staged decoupling schedule**: Electrostatic point charges must be turned off linearly before scaling Lennard-Jones van der Waals terms using soft-core potentials. This staged protocol prevents unshielded point charges from overlapping with surrounding solvent or protein atoms (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).
- **Cycle closures and diagnostics**: Discrepancies between forward and backward work distributions (hysteresis) or failure of thermodynamic cycle closure ($\sum \Delta\Delta G_{\text{cycle}} \ne 0$) diagnose hidden barriers or insufficient sampling (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf). Dedicated details are covered in [[alchemical-free-energy-methods]].

## Practical domain expansions and machine learning

Recent studies have expanded the domain of applicability of alchemical calculations (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf; Integrating Machine Learning into Free Energy Perturbation Workflows.pdf):
- **Homology models**: Across 29 homology models of 14 protein targets, FEP results on models with 30%–50% template sequence identity matched crystal structure accuracy after induced-fit MD relaxation (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).
- **Fragment growing and covalent inhibitors**: Evaluated on 90 fragments across eight targets to predict hit elaboration, and utilized to compute the non-covalent binding step ($K_i$) of covalent inhibitors prior to bond formation ($k_{\text{inact}}$) (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).
- **Active learning (AL) acceleration**: AL coupled with surrogate graph neural networks or Gaussian processes screens multi-million-compound libraries up to 20-fold faster by computing FEP only on high-utility acquisition candidates (source: Integrating Machine Learning into Free Energy Perturbation Workflows.pdf).
- **Machine learning potentials (MLPs)**: Neural network potentials (AIM-Net2, MACE) deployed via hybrid ML/MM (DP/MM) simulations replace classical fixed-charge force fields, providing quantum chemical torsional and polarization accuracy within the alchemical cycle (source: Integrating Machine Learning into Free Energy Perturbation Workflows.pdf).
- **Endpoint methods**: In high-throughput virtual screening, [[mm-pbsa-and-mm-gbsa]] serves as a critical intermediate filter, evaluating single-trajectory explicit MD snapshots to filter out false positives and boost early enrichment factors (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf).

## What limits accuracy

1. **Structural input.** Historically, high-resolution co-crystal structures (<2.5 Å) were considered mandatory (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). However, recent large-scale benchmarks demonstrate that homology models (30%–50% sequence identity) and AlphaFold/Chai-1 holo models achieve predictive accuracy comparable to experimental crystal structures ($\text{RMSE} \sim 1.0\text{--}1.5\text{ kcal/mol}$) once sidechains and pocket loops are relaxed through initial equilibrium MD simulations (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).
2. **Force field and charges.** Charge-changing perturbations and charged ligands are especially sensitive (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
3. **Sampling and ergodicity.** Slow conformational modes can be missed; poor ergodicity can arise in protein, ligand, or solvent (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulation for all.md).
4. **Standard-state and restraint corrections.** Ignoring them makes absolute values incomparable with experiment.
5. **Water and entropy.** Both are difficult and can dominate (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
6. **Prospective validation.** Most studies are retrospective; blind challenges and common benchmarks are still limited (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

## Practical protocol

1. Start from a high-quality structure or well-relaxed homology/AlphaFold model; generate poses by docking only if no co-crystal exists, then cluster and refine with MD equilibration (source: Molecular dynamics simulation for all.md; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).
2. Build both end-state topologies; choose single or dual topology and a soft-core scheme (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf).
3. Use a validated force field, consistent charges, and explicit solvent, optionally applying [[hydrogen-mass-repartitioning]] to permit 4.0 fs timesteps (source: Molecular dynamics simulation for all.md).
4. Run multiple λ windows with adequate equilibration; use FEP/REST or Hamiltonian replica exchange when barriers are high (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
5. Apply restraints for absolute calculations and compute the analytical restraint correction (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
6. Repeat with independent seeds; report mean, standard deviation, and convergence.
7. Validate against [[experimental-binding-assays]] (ITC or SPR data) on a congeneric series before prospective deployment (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

## Related pages
- [[alchemical-free-energy-methods]]
- [[mm-pbsa-and-mm-gbsa]]
- [[binding-thermodynamics]]
- [[binding-kinetics]]
- [[molecular-dynamics-simulation]]
- [[hydrogen-mass-repartitioning]]
- [[enhanced-sampling]]
- [[force-fields]]
- [[water-and-solvation]]
- [[scoring-functions]]
- [[experimental-binding-assays]]
- [[decherchi-2020-binding-thermodynamics-kinetics]]
- [[genheden-2015-mm-pbsa-gbsa]]
- [[muegge-2023-alchemical-free-energy-advances]]
- [[york-2023-modern-alchemical-free-energy]]
- [[adediwura-2024-free-energy-kinetics-drug-discovery]]
- [[metadynamics]]
- [[umbrella-sampling]]
