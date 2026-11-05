# Simulation and Benchmark Best Practices

**Summary**: A consolidated best-practices guide synthesizing the LiveCoMS guidelines for foundational molecular dynamics simulation validation and the Hahn et al. standards for constructing, preparing, and evaluating protein–ligand binding affinity benchmarks.
**Sources**: `Best Practices for Foundations in Molecular Simulations [Article.pdf`, `Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf`, `mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf`, `Assessment of molecular dynamics time series.pdf`
**Last updated**: 2026-11-02

---

## Foundational guidelines for molecular dynamics simulations

Biomolecular simulations produce physically meaningful results only when statistical mechanical ensembles, numerical integrators, boundary conditions, and thermodynamic regulators are implemented rigorously (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). The Living Journal of Computational Molecular Science (LiveCoMS) established foundational standards to prevent common simulation artifacts and methodological failures (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

### 1. Integrator Selection and Timestep Boundaries
- **Timestep Limit**: The integration timestep $\Delta t$ must never exceed one-tenth of the period of the fastest vibrational mode ($\Delta t \le \tau_{\text{fastest}} / 10$) (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
  - Unconstrained all-atom simulations: $\Delta t \le 1.0\text{ fs}$ (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
  - Constrained all-atom simulations (SHAKE/LINCS on hydrogen bonds): $\Delta t = 2.0\text{ fs}$ (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
  - Hydrogen Mass Repartitioning ([[hydrogen-mass-repartitioning]]) or virtual interaction sites (VSITE): enables $\Delta t = 4.0\text{ fs}$ without altering configuration-space thermodynamic ensembles (source: Best Practices for Foundations in Molecular Simulations [Article.pdf; source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
- **Microcanonical ($NVE$) Drift Check**: Any chosen integration timestep and cutoff scheme must be validated by running a short $NVE$ trajectory; the total energy should exhibit zero systematic upward or downward drift (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

### 2. Thermostat Audit and Known Artifacts
- **Algorithms to Avoid**: Simple velocity rescaling and the **Berendsen thermostat** must not be used in production simulations (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). They violate the equipartition theorem, fail to generate a true canonical ($NVT$) ensemble, and trigger the notorious **"flying ice cube"** artifact, where kinetic energy drains from high-frequency vibrational modes into rigid center-of-mass translation (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Recommended Thermostats**: The Bussi stochastic velocity-rescaling thermostat (`v-rescale`) or the **Nosé–Hoover** chain thermostat (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). Bussi combines velocity rescaling with a Wiener stochastic noise process, strictly enforcing the canonical probability distribution while maintaining momentum conservation and physical dynamics (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Temperature Coupling Groups**: Small molecules and bound ligands must never be thermostat-coupled independently or lumped with bulk solvent; they must be merged with the protein (`tc-grps = Protein_Ligand Water_and_ions`) to prevent spurious thermal gradients (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

### 3. Barostats and Isobaric Ensembles
- **Equilibration vs. Production**: The Berendsen barostat is acceptable strictly during early equilibration due to its numerical robustness, but it suppresses volume fluctuations and must never be used for production MD or property calculations (compressibility, free energy) (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Rigorous $NPT$ Barostats**: Production simulations require the **Parrinello–Rahman**, **Andersen**, or **Monte Carlo** barostats, which generate the true isothermal–isobaric ensemble (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

### 4. Periodic Boundary Conditions and Nonbonded Cutoffs
- **Minimum Image Convention**: To prevent a macromolecule from interacting with its own periodic mirror images, the shortest dimension of the simulation cell must be at least twice the nonbonded cutoff distance ($L_{\text{box}} \ge 2 r_{\text{cut}}$) (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Cell Geometry**: Rhombic dodecahedral and triclinic boxes should be preferred over cubic boxes for globular systems; they reduce explicit solvent volume by $\sim 23\%$ for an identical solute-to-wall distance, reducing compute cost substantially (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
- **Long-Range Electrostatics**: Particle Mesh Ewald (PME) electrostatics with a grid spacing $\le 0.16\text{ nm}$ and fourth-order spline interpolation is mandatory (source: Best Practices for Foundations in Molecular Simulations [Article.pdf; source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- **Dispersion Corrections**: Van der Waals interactions require force-switching (CHARMM) or continuous cutoffs with analytical dispersion corrections (AMBER) to prevent density and pressure artifacts (source: Best Practices for Foundations in Molecular Simulations [Article.pdf; source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

### 5. Replicate Sampling and Convergence Verification
- **Multiple Independent Trajectories**: Sampling convergence cannot be established from a single long simulation (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). Running $3\text{ to }5$ independent replicate simulations initialized with different random Maxwell–Boltzmann velocity seeds is mandatory to detect kinetic trapping and calculate accurate statistical error bars (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Equilibration Trimming ("Burn-in")**: Properties such as solvent density, hydrogen-bond networks, and pocket conformations require distinct relaxation times; initial non-equilibrium phases must be discarded prior to ensemble averaging (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

---

## Best practices for binding affinity benchmark construction

Evaluating computational drug discovery tools—such as relative binding free energy (RBFE), absolute binding free energy (ABFE), or machine learning scoring functions—requires standardized benchmarks adherent to the LiveCoMS guidelines established by Hahn et al. (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

### 1. Experimental Assay Curation Criteria
- **Direct Biophysical Assays**: Prioritize direct thermodynamic measurements—Isothermal Titration Calorimetry (ITC) and Surface Plasmon Resonance (SPR)—which determine true equilibrium dissociation constants ($K_d$) independent of enzymatic substrate competition (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Handling $\text{IC}_{50}$ Data**: Avoid functional or cellular $\text{IC}_{50}$ values where possible (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Conversion to $K_i$ via Cheng–Prusoff is valid only when the enzyme concentration is far below the inhibition constant ($[E] \ll K_i$) and substrate saturation kinetics are rigorously defined (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Dynamic Range**: A benchmark set must exhibit an experimental affinity span of at least $2\text{ to }3$ orders of magnitude ($\Delta\Delta G \ge 2.7\text{--}4.1\text{ kcal/mol}$) to provide statistical power for correlation metrics (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Experimental Noise Floor**: Replicate experimental measurements exhibit an inherent uncertainty of $\sim 0.3\text{ to }0.5\text{ kcal/mol}$ (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Computational predictions achieving errors within this window have reached the experimental resolution limit (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

### 2. Structural Quality and System Preparation
- **Crystallographic Metrics**: High-resolution X-ray structures ($\le 2.0\text{--}2.5\text{ \AA}$) with low free $R$-factors ($R_{\text{free}} < 0.25$) and high Real-Space Correlation Coefficients (RSCC) are required (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Electron density maps ($2F_o - F_c$) must confirm unambiguous placement of binding-pocket sidechains and ligand functional groups (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Crystal Packing Inspection**: Verify that active-site loops are not stabilized or distorted by crystal lattice contacts with neighboring symmetry mates (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Protonation and Tautomer Modeling**: Amino acid ionization states (His, Asp, Glu, Lys) and ligand tautomers must be modeled at the exact pH of the experimental assay (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Bridging Waters and Metals**: Retain functionally conserved water molecules and catalytic metals (e.g., $\text{Zn}^{2+}$, $\text{Mg}^{2+}$) with verified coordination geometries (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

### 3. Congeneric Series Design for Alchemical Free Energy (RBFE)
- **Closed-Loop Perturbation Graphs**: Networks must contain redundant, closed thermodynamic cycles (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Cycle closure hysteresis ($\sum_{\text{cycle}} \Delta\Delta G = 0$) provides an automated, intrinsic metric of simulation convergence (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Maximum Common Substructure (MCS)**: Transformations between ligand pairs should minimize the number of alchemically perturbed non-hydrogen atoms, avoiding large scaffold hops, ring openings/closings, or net charge changes that introduce phase-space sampling barriers (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

### 4. Statistical Evaluation Best Practices
- **Misleading Nature of $R^2$**: Do not rely exclusively on the coefficient of determination ($R^2$); it artificially inflates when testing broad dynamic ranges and obscures large localized errors (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Mandatory Reporting Metrics**: Report root-mean-square error (RMSE), mean absolute error (MAE / MUE), and rank-order correlation metrics—specifically Kendall's $\tau$ (which evaluates pairwise concordant vs. discordant ranking) and Spearman's $\rho$ (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Bootstrapping**: All performance metrics must include non-parametric bootstrapped $95\%$ confidence intervals (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

---

## Trajectory feature curation for machine learning models

When integrating molecular dynamics trajectories into machine learning scoring functions:
- **Simulation Noise Hazard**: Training ML models directly on uncurated trajectory frames degrades performance as trajectory length increases, because thermal fluctuations introduce high-dimensional noise (source: Assessment of molecular dynamics time series.pdf).
- **Multi-Stage Feature Selection**: Time-series descriptors (`ts_descriptors`) extracted from trajectories must undergo systematic pruning:
  1. $p$-value filtering against experimental binding affinity (source: Assessment of molecular dynamics time series.pdf).
  2. Multi-collinearity removal (pruning pairs with Pearson correlation coefficient $\text{PCC} \ge 0.8$) (source: Assessment of molecular dynamics time series.pdf).
  3. Mutual information maximization to select the single most informative descriptor per biophysical property (source: Assessment of molecular dynamics time series.pdf).
- **Generalization Regimes**: Curated MD descriptors specifically enhance out-of-distribution generalizability on challenging **Target Splits** (unseen protein families) and **Scaffold Splits** (unseen ligand chemotypes) by capturing physical interaction persistence over time (source: Assessment of molecular dynamics time series.pdf).

## Related pages
- [[braun-2019-foundations-molecular-simulation]]
- [[hahn-2022-binding-affinity-benchmarks]]
- [[gromacs-protein-ligand-protocol]]
- [[molecular-dynamics-simulation]]
- [[benchmarking-datasets]]
- [[binding-free-energy]]
- [[alchemical-free-energy-methods]]
- [[open-force-field-and-ligand-parameterization]]
- [[poziemski-2026-md-time-series-descriptors]]
- [[machine-learning-scoring-functions]]

