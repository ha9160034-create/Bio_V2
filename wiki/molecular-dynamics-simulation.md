# Molecular Dynamics Simulation

**Summary**: Molecular dynamics (MD) integrates Newton's equations of motion for every atom in a solvated biomolecular system using a molecular-mechanics force field, producing a trajectory that describes structure, dynamics, and thermodynamics at atomic resolution.
**Sources**: Molecular dynamics simulation for all.md; Molecular dynamics simulations and drug discovery.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf; Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Enhanced Sampling Techniques in Molecular Dynamics Simulations of Biological Systems.pdf; The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf; Best Practices for Foundations in Molecular Simulations [Article.pdf; gromacs_protein_ligand_complex_tutorial_lemkul.md; mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf; Assessment of molecular dynamics time series.pdf
**Last updated**: 2026-11-02

---

## The basic algorithm

Given the positions of all atoms in a system — protein, ligand, water, ions, and optionally a lipid bilayer — the force on each atom is computed from the force field. Newton's laws then update each position and velocity, the time is advanced by a small step, and the process repeats (source: Molecular dynamics simulation for all.md).

The result is a trajectory: a three-dimensional movie of the atomic configuration at every point in the simulated interval (source: Molecular dynamics simulation for all.md).

MD is powerful because it captures the position and motion of every atom at every time point, and because the conditions are precisely controlled and can be perturbed: initial conformation, bound ligands, mutations, post-translational modifications, solvent composition, protonation states, temperature, and membrane voltage (source: Molecular dynamics simulation for all.md).

## Time step and cost

To keep numerical integration stable, the time step is typically a few femtoseconds (10⁻¹⁵ s). Most biochemical events of interest occur on nanoseconds to microseconds or longer, so a simulation requires millions to billions of steps. Each step evaluates millions of interatomic interactions, which makes MD computationally demanding (source: Molecular dynamics simulation for all.md).

The common practical choice is **1 fs**, or **2 fs when bonds involving hydrogen are constrained** (source: Molecular dynamics simulation for all.md). The Vina/MD literature likewise notes that a 1 fs step means a millisecond event requires 10¹² integration steps (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

**Hydrogen mass repartitioning (HMR)** transfers atomic mass from heavy atoms (C, N, O) to attached hydrogens, raising hydrogen mass to ~3.0 amu (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf):
- Slows high-frequency bond-stretching modes without altering the potential energy surface;
- Expands the stable numerical integration timestep to **4.0 fs**, granting a **two-fold acceleration in sampling throughput** while preserving thermodynamic ensembles and order parameters (Hopkins et al., 2015; Balusek et al., 2019) (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). See [[hydrogen-mass-repartitioning]].

Historical scale puts the cost in perspective: in 2011, a 1 µs simulation of about 25,000 atoms on 24 processors took several months (source: Molecular dynamics simulations and drug discovery.md). By 2018, one or two GPUs could outperform what previously required most of a supercomputer, and a 1 µs simulation of a ~50,000-atom system could be run in a few days (source: Molecular dynamics simulation for all.md). The two statements are not contradictory; they reflect a decade of hardware and software acceleration.

## What MD can answer

MD is used for many question types (source: Molecular dynamics simulation for all.md):

1. **Flexibility and mobility** — quantify how much regions of a structure move at equilibrium, and how water and ions behave.
2. **Structure validation and refinement** — relax crystal-packing or membrane-absence artifacts; test whether a docked pose is stable; refine cryo-EM models with ambiguous ligand density.
3. **Model building from experimental data** — MD-based simulated annealing for crystallographic refinement, flexible fitting into cryo-EM density, and ensemble recovery from NMR restraints.
4. **Perturbation response** — remove, replace, or add a ligand; mutate a residue; phosphorylate; change protonation; apply voltage or mechanical strain; change salt or lipid composition.
5. **Mechanism** — observe ligand binding, conformational change, folding, allostery, and transport in action.
6. **Thermodynamics and kinetics** — estimate binding free energies and rates, usually with enhanced sampling or alchemical methods.

## Drug-discovery applications

### Cryptic and allosteric site discovery

MD revealed the cryptic trench on HIV integrase that led to raltegravir, and FTMAP on MD conformations of β-adrenergic receptors identified five potential allosteric sites, some absent from crystal structures (source: Molecular dynamics simulations and drug discovery.md). Furthermore, [[mixed-solvent-molecular-dynamics]] (MixMD, SILCS, SWISH) sweeps macromolecular surfaces with small organic probe mixtures (e.g., 5% v/v isopropanol, benzene, or phenol) to compete against water, inducing cryptic pocket opening and generating grid free energy (GFE) hotspot maps for allosteric drug design (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). See [[cryptic-and-allosteric-pockets]] and [[mixed-solvent-molecular-dynamics]].

### Ensemble docking and the relaxed complex scheme

Instead of docking into one static structure, the relaxed complex scheme docks each ligand into multiple MD-derived receptor conformations, producing a spectrum of scores per ligand. Ligands are ranked by a statistic of that spectrum, such as the mean score. RCS has identified inhibitors of FKBP, HIV integrase, *Trypanosoma brucei* RNA editing ligase 1, *T. brucei* GalE, *T. brucei* FPPS, and *M. tuberculosis* dTDP-6-deoxy-L-lyxo-4-hexulose synthase (source: Molecular dynamics simulations and drug discovery.md).

### Binding free energy

Alchemical free-energy methods — thermodynamic integration, free-energy perturbation, and single-step perturbation — are built on MD trajectories. Relative free energies for congeneric series and absolute free energies via double decoupling are the most accurate computational affinity estimates available, at high computational cost (source: Molecular dynamics simulations and drug discovery.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[binding-free-energy]].

### Binding kinetics

Unbiased multi-microsecond MD can observe spontaneous binding and unbinding and yield k_on and k_off. Markov state models built from many short trajectories extend the accessible timescale, and weighted ensemble methods reach residence times far beyond plain MD (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[binding-kinetics]].

### Pose validation and recognition tracking

A docked pose that remains stable in MD is more likely to be correct than one that rapidly distorts. This has been used to resolve ligand placement in cryo-EM structures with ambiguous density (source: Molecular dynamics simulation for all.md). Lemkul's standardized setup in [[gromacs-protein-ligand-protocol]] demonstrates that positional restraints during $NVT$/$NPT$ equilibration, combined with merged temperature coupling groups (`Protein_JZ4`), maintain native complex integrity while solvent relaxes (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). Furthermore, bridging docking with MD enables active tracking of the binding pathway; Supervised Molecular Dynamics (SuMD) evaluates short classical MD windows (e.g., 600 ps) to monitor ligand approach and binding pocket entry without unphysical bias forces (source: Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf). In post-docking virtual screening, short explicit MD trajectories (0.5–5 ns) provide conformational ensembles for endpoint [[mm-pbsa-and-mm-gbsa]] free-energy calculations, filtering out unstable docking artifacts and boosting screening enrichment (source: The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf). For high-barrier conformational changes, generalized simulated annealing and non-equilibrium enhanced sampling expand phase-space exploration (source: Enhanced Sampling Techniques in Molecular Dynamics Simulations of Biological Systems.pdf). See [[enhanced-sampling]] and [[gromacs-protein-ligand-protocol]].

### Machine learning trajectory descriptors

Extracting dynamic time-series features (`ts_descriptors`) from MD trajectories overcomes the fundamental limitation of static crystal structures in affinity prediction (source: Assessment of molecular dynamics time series.pdf). On a benchmark of 2,752 complexes, Poziemski et al. demonstrated that raw, unpruned MD frames introduce simulation noise that degrades ML scoring, whereas aggressive multi-stage feature trimming ($p$-value filtering, correlation pruning at $\text{PCC} \ge 0.8$, and mutual information maximization) captures essential contact persistence and energy fluctuations that significantly improve model generalization on challenging Target Splits and Scaffold Splits (source: Assessment of molecular dynamics time series.pdf). See [[poziemski-2026-md-time-series-descriptors]] and [[machine-learning-scoring-functions]].

## Force fields and solvent

The forces come from a molecular-mechanics force field fit to quantum-mechanical calculations and experimental data. Typical terms are Coulombic electrostatics, spring-like bond stretching, angle bending, sinusoidal torsions, and Lennard-Jones van der Waals interactions (source: Molecular dynamics simulation for all.md; Molecular dynamics simulations and drug discovery.md).

Classical MD cannot break or form covalent bonds; QM/MM is used when chemistry occurs in a small region (source: Molecular dynamics simulation for all.md). See [[force-fields]].

Water is usually explicit. The choice of water model, ion parameters, and membrane composition materially affects results, and the force field remains approximate despite substantial improvement over the past decade (source: Molecular dynamics simulation for all.md).

## Limitations

- **Force-field accuracy.** Force fields are imperfect approximations; uncertainty must be considered when interpreting results. Polarization is usually ignored, and a generally accepted polarizable force field has not emerged (source: Molecular dynamics simulations and drug discovery.md; Molecular dynamics simulation for all.md).
- **Sampling.** Many processes occur on timescales far longer than routine simulations. Insufficient sampling can miss relevant conformations entirely (source: Molecular dynamics simulations and drug discovery.md).
- **Quantum effects.** Transition-metal chemistry, bond breaking, and electronic polarization require QM/MM or specialized methods (source: Molecular dynamics simulations and drug discovery.md).
- **Interpretation.** Producing a trajectory is easy; deciding which questions it can answer and extracting biological insight from a mass of jiggling atoms is hard (source: Molecular dynamics simulation for all.md).

## Hardware and software

Popular packages include GROMACS, NAMD, AMBER, CHARMM, Desmond, and OpenMM. GPUs and specialized hardware such as Anton have pushed simulations to millisecond timescales for selected systems (source: Molecular dynamics simulation for all.md; Molecular dynamics simulations and drug discovery.md).

## Practical guidance

- Start from a high-quality structure. While historical protocols strictly required experimental crystal structures better than 2.5 Å (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md), contemporary benchmarks establish that homology models (30%–50% identity) and AlphaFold holo models achieve quantitative accuracy comparable to experimental co-crystals after initial solvent equilibration and MD relaxation (source: Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf; Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).
- Follow foundational LiveCoMS simulation checklists: enforce minimum image conventions ($L_{\text{box}} \ge 2 r_{\text{cut}}$), utilize rhombic dodecahedral boxes to reduce solvent water volume by $\sim 23\%$ over cubic boxes, avoid the Berendsen thermostat in production runs to prevent the "flying ice cube" artifact, and select Bussi stochastic velocity rescaling (`v-rescale`) or Nosé–Hoover with Parrinello–Rahman or Monte Carlo barostats (source: Best Practices for Foundations in Molecular Simulations [Article.pdf; mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf). See [[simulation-and-benchmark-best-practices]].
- Run several independent replicates of both perturbed and unperturbed systems with different initial velocity seeds to identify consistent differences and avoid kinetic trapping (source: Molecular dynamics simulation for all.md; Best Practices for Foundations in Molecular Simulations [Article.pdf).
- Check convergence and replicate variance; a single trajectory can differ from another by as much as the full range of interest (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- Prefer methods with few user-tuned choices; methods whose outcome depends heavily on human choices are not mature enough for prospective quantitative use (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Related pages
- [[force-fields]]
- [[enhanced-sampling]]
- [[hydrogen-mass-repartitioning]]
- [[gromacs-protein-ligand-protocol]]
- [[mixed-solvent-molecular-dynamics]]
- [[simulation-and-benchmark-best-practices]]
- [[binding-free-energy]]
- [[mm-pbsa-and-mm-gbsa]]
- [[alchemical-free-energy-methods]]
- [[binding-thermodynamics]]
- [[binding-kinetics]]
- [[water-and-solvation]]
- [[cryptic-and-allosteric-pockets]]
- [[benchmark-drug-targets]]
- [[poziemski-2026-md-time-series-descriptors]]
- [[braun-2019-foundations-molecular-simulation]]
- [[lemkul-gromacs-protein-ligand-tutorial]]
- [[lemkul-2018-gromacs-tutorials]]
- [[hollingsworth-2018-md-for-all]]
- [[durrant-2011-md-drug-discovery]]
- [[decherchi-2020-binding-thermodynamics-kinetics]]
- [[bernardi-2015-enhanced-sampling-md]]
- [[salmaso-2018-bridging-docking-md]]
