# Hollingsworth & Dror (2018) — MD Simulation for All

**Summary**: A review explaining what molecular dynamics simulations are, how they are run, what they can reveal about biomolecular function and drug binding, and why recent hardware and software advances have brought the method to a much wider audience.
**Sources**: Molecular dynamics simulation for all.md; Molecular dynamics simulations and drug discovery.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md
**Last updated**: 2026-11-01

---

## Citation and scope

Sarah A. Hollingsworth and Ron O. Dror, "Molecular Dynamics Simulation for All," *Neuron* 99(6):1129–1143 (2018), doi:10.1016/j.neuron.2018.08.011. PMC6209097 (source: Molecular dynamics simulation for all.md).

The review is written for experimental biologists and clinicians rather than simulation specialists. It explains MD conceptually, describes practical choices, surveys applications across molecular biology and drug discovery, and is candid about limitations. It is the wiki's main general-purpose MD reference.

## The central idea

Given the positions of all atoms in a biomolecular system, one computes the force on each atom from the force field and uses Newton's laws to update positions and velocities; the cycle repeats at a small time step. The output is a trajectory — in effect a three-dimensional movie of the atomic configuration over the simulated interval (source: Molecular dynamics simulation for all.md).

MD is powerful for three reasons (source: Molecular dynamics simulation for all.md):

1. it captures the position and motion of every atom at every time point, which no experiment does routinely;
2. the conditions are precisely known and controlled — starting conformation, bound ligands, mutations, post-translational modifications, solvent, protonation states, temperature, membrane voltage;
3. the system can be perturbed at will and the response observed.

The review's analogy is a molecular biologist watching a protein work rather than examining a single still photograph: function and interactions depend on motion, not only on structure (source: Molecular dynamics simulation for all.md).

## Time step and computational cost

The time step must be short for numerical stability — typically a few femtoseconds (10⁻¹⁵ s). Biochemical events of interest occur on nanoseconds to microseconds or longer, so a simulation requires millions to billions of steps, each evaluating millions of interatomic interactions (source: Molecular dynamics simulation for all.md).

The practical choices are **1 fs**, or **2 fs when bonds involving hydrogen are constrained** (source: Molecular dynamics simulation for all.md).

By 2018, a 50,000-atom solvated protein could be simulated for 1 µs in a few days on one GPU, while specialized hardware could increase speed by at least an order of magnitude at much higher cost (source: Molecular dynamics simulation for all.md). This is a dramatic change from the 2011 figure of several months on 24 processors for a similar system (source: Molecular dynamics simulations and drug discovery.md); the two are consistent once hardware and software advances are accounted for.

## Force field

Forces come from a molecular-mechanics force field fitted to quantum-mechanical calculations and, typically, to experimental measurements. A typical force field contains Coulombic electrostatics, spring-like bond-length terms, angle terms, torsional terms, and Lennard-Jones van der Waals interactions. The functional form is approximate, and force fields have improved substantially over the past decade but remain imperfect; the uncertainty should be considered when interpreting results (source: Molecular dynamics simulation for all.md).

Classical MD cannot break or form covalent bonds. **QM/MM** treats a small chemically active region quantum mechanically and the remainder with molecular mechanics (source: Molecular dynamics simulation for all.md). See [[force-fields]].

The review highlights three force-field families (source: Molecular dynamics simulation for all.md):

- **CHARMM36m** and the CHARMM General Force Field (CGenFF) — extensively optimized for proteins, lipids, and drug-like ligands;
- **A99SB-disp** — particularly good for intrinsically disordered proteins;
- **OPLS3** — perhaps the most extensively optimized ligand parameters, although protein parameters are less thoroughly validated.

Water is usually explicit, and the water model and ion parameters must be consistent with the force field (source: Molecular dynamics simulation for all.md). See [[water-and-solvation]].

## Software and hardware

Popular packages include **GROMACS, NAMD, AMBER, CHARMM, Desmond, and OpenMM** (source: Molecular dynamics simulation for all.md). GPUs have been the decisive hardware change: one or two inexpensive GPUs can outperform what previously required most of a supercomputer. Specialized machines such as Anton reach millisecond timescales for selected systems (source: Molecular dynamics simulation for all.md; Molecular dynamics simulations and drug discovery.md).

## What MD can answer

The review organizes applications into six groups (source: Molecular dynamics simulation for all.md):

1. **Flexibility and mobility** — how much regions of a structure move at equilibrium, and how water and ions behave.
2. **Structure validation and refinement** — relax crystal-packing or membrane-absence artifacts, test whether a docked pose is stable, and refine cryo-EM models with ambiguous ligand density.
3. **Model building from experimental data** — MD-based simulated annealing for crystallographic refinement, flexible fitting into cryo-EM density, and ensemble recovery from NMR restraints.
4. **Perturbation response** — remove, replace, or add a ligand; mutate a residue; phosphorylate; change protonation; apply voltage or mechanical strain; change salt or lipid composition.
5. **Mechanism** — observe ligand binding, conformational change, folding, allostery, and transport.
6. **Thermodynamics and kinetics** — estimate binding free energies and rates, usually with enhanced sampling or alchemical methods.

## Applications in drug discovery

### Cryptic and allosteric sites

MD can reveal pockets absent from crystal structures. The canonical example is HIV integrase: simulations exposed a trench that later crystallography showed to be the binding site of inhibitors, work that contributed to raltegravir (source: Molecular dynamics simulations and drug discovery.md). The review discusses analogous uses of MD to find allosteric sites and to explain allosteric mechanisms (source: Molecular dynamics simulation for all.md). See [[cryptic-and-allosteric-pockets]].

### Docking and pose validation

A docked pose that remains stable during MD is more likely to be correct than one that rapidly distorts. MD has been used to resolve ligand placement in cryo-EM structures with ambiguous density and to relax structures before docking (source: Molecular dynamics simulation for all.md). See [[protein-ligand-docking]].

### Binding free energy and kinetics

Alchemical free-energy methods are built on MD trajectories; they give the most accurate computational affinity estimates available, at high cost. Unbiased long simulations and enhanced sampling can also estimate binding and unbinding rates (source: Molecular dynamics simulation for all.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[binding-free-energy]] and [[binding-kinetics]].

### Allostery

Simulations of the M2 muscarinic acetylcholine receptor showed how allosteric modulators change the conformation of the orthosteric site. The orthosteric ligand N-methyl scopolamine favors an enlarged allosteric site; the positive modulator alcuronium requires the larger site, whereas the negative modulator does not. A designed modulator predicted to need a larger pocket had fourfold less negative cooperativity despite higher affinity — an experimental validation of a simulation-derived mechanism (source: Molecular dynamics simulation for all.md).

### Membrane proteins and channels

MD is used to study transporters, channels, and GPCRs in a lipid bilayer, including ion permeation, lipid modulation, and the effect of mutations (source: Molecular dynamics simulation for all.md).

## Limitations

- **Force-field accuracy.** Force fields are approximate, polarization is usually ignored, and no generally accepted polarizable force field exists (source: Molecular dynamics simulation for all.md; Molecular dynamics simulations and drug discovery.md).
- **Sampling.** Many processes occur on timescales longer than routine simulations; insufficient sampling can miss relevant conformations entirely (source: Molecular dynamics simulation for all.md).
- **Quantum effects.** Bond breaking, transition-metal chemistry, and electronic polarization require QM/MM or specialized methods (source: Molecular dynamics simulation for all.md; Molecular dynamics simulations and drug discovery.md).
- **Interpretation.** Producing a trajectory is straightforward; deciding which questions it can answer and extracting biological insight from a mass of jiggling atoms is hard (source: Molecular dynamics simulation for all.md).
- **Force-field dependence of mechanism.** A simulated mechanism can be an artifact of the parameter set; independent replicates and experimental tests are essential.

## Practical guidance

The review's advice for running and interpreting simulations (source: Molecular dynamics simulation for all.md):

1. Choose hardware according to the timescale required; GPUs are usually the best price-performance choice.
2. Prepare the system carefully: build missing atoms and hydrogens, add water, ions, and lipids, and assign protonation states.
3. Run **multiple independent simulations** under each condition, typically starting from the same coordinates with randomized initial velocities, to distinguish real differences from trajectory-to-trajectory variation.
4. Track a well-defined quantity rather than trying to absorb the entire trajectory.
5. Compare with experimental data — mutagenesis, spectroscopy, NMR relaxation, or crystallographic B-factors.
6. Consider enhanced sampling or coarse-graining when the process of interest is beyond accessible timescales, accepting the associated loss of accuracy.

## Why it matters for this wiki

This is the main source for [[molecular-dynamics-simulation]] and [[force-fields]], and a supporting source for [[cryptic-and-allosteric-pockets]], [[binding-free-energy]], [[binding-kinetics]], and [[water-and-solvation]]. It provides the wiki's operational details — time step, force-field families, software, hardware, and the six application classes — that make the MD pages concrete.

## Related pages
- [[molecular-dynamics-simulation]]
- [[force-fields]]
- [[enhanced-sampling]]
- [[binding-free-energy]]
- [[binding-kinetics]]
- [[cryptic-and-allosteric-pockets]]
- [[water-and-solvation]]
- [[durrant-2011-md-drug-discovery]]
- [[decherchi-2020-binding-thermodynamics-kinetics]]
