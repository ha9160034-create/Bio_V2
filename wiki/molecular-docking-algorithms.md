# Molecular Docking Algorithms

**Summary**: Docking algorithms search the conformational space of a ligand (and sometimes a receptor) using systematic, stochastic, simulation-based, or machine-learning strategies, trading completeness against computational cost.
**Sources**: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Protein–Ligand Docking in the Machine-Learning Era.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Basic docking — Autodock Vina 1.2.0 documentation.md
**Last updated**: 2026-11-01

---

## The search problem

The theoretical search space for protein–ligand binding includes all conformations of the unbound protein and ligand, all ligand orientations and conformations within each protein state, and all protein conformations paired with all ligand conformations. Exhaustive exploration is impossible with current computers (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

A ligand contributes six rigid-body degrees of freedom (three translations, three rotations) plus one torsional degree of freedom per rotatable bond. Proteins add backbone and side-chain dihedrals, and in solution they exist as an ensemble of interconverting conformations (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). The search algorithm must therefore balance **speed** and **coverage** of the relevant conformational space (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Rigid-body docking

The simplest algorithms treat both partners as rigid bodies and explore only the six rotational and translational degrees of freedom of the ligand. Examples include ZDOCK, MDock, and older versions of DOCK (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

“Rigid-body” is not identical to “rigid ligand”: ligand flexibility is often handled by pre-computing an ensemble of ligand conformations and docking each conformer rigidly (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Flexible-ligand docking

Flexible-ligand algorithms keep the receptor rigid but generate ligand conformations during the search. Three major strategies exist (source: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md):

### Systematic search and incremental construction

Systematic methods enumerate torsions, which risks exponential growth of the search space. Grid-based storage of physicochemical properties mitigates combinatorial explosion. The ligand is divided into fragments that are docked sequentially and then assembled — the **anchor-and-grow** or incremental construction (IC) strategy (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

- DOCK, the first and most widely used docking program, implements incremental construction; the current version is DOCK 6 (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- FlexX and Glide also use incremental construction (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- LUDI builds ligands from fragments placed in the binding cavity (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Stochastic search

Stochastic methods introduce randomness to escape local minima (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md):

- **Monte Carlo (MC) / simulated annealing**: random changes in position, orientation, or conformation are accepted or rejected by the Metropolis criterion. The simulation starts at high temperature and cools gradually (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Genetic algorithms (GA)**: a population of ligand states is encoded as chromosomes; crossover and mutation generate new generations, and the lowest-energy individuals survive. AutoDock uses a GA, and AutoDockFR combines a GA with a Solis–Wets local search (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### Molecular-dynamics-based search

MD can be used for docking sampling, but a plain MD trajectory tends to become trapped in local minima because a move to a higher-energy state is disfavoured. A common compromise is to identify candidate conformations first and then use MD for local optimization; raised temperature and smoothed potential surfaces are also used (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## Receptor flexibility

Most docking programs treat the protein as rigid because of its size. Several strategies approximate receptor flexibility (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Protein–Ligand Docking in the Machine-Learning Era.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- **Ensemble docking**: dock the ligand into multiple receptor conformations from X-ray, NMR, or MD, and select the best-scoring complex. MDock treats the protein conformational number as an additional variable, making ensemble docking roughly as fast as single-structure docking. Excessive numbers of conformers increase false positives and cost linearly.
- **Soft docking**: allow limited steric overlap between ligand and protein atoms, or optimize selected side-chain rotamers. AutoDockFR samples explicitly specified side chains with the AutoDock4 force field; it improves pose prediction but can increase false positives in virtual screening.
- **Enlarged binding pockets**: ADAM shifts van der Waals curves to widen the cavity; SCARE enlarges the pocket by alanine mutation.
- **Partial side-chain flexibility**: rotamer libraries or flexible side chains sampled simultaneously with the ligand (AutoDock, ICM, GalaxyDock).
- **Energy minimization**: MC or gradient-descent minimization of the complex, corresponding conceptually to induced fit.

## Machine-learning pose prediction

Deep-learning methods now predict poses directly, often in a single shot (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md):

- **EquiBind** uses SE(3)-equivariant graph neural networks and treats docking as regression, predicting pocket keypoints on protein and ligand and aligning them.
- **TANKBind** adds trigonometry constraints and segments the protein into functional blocks to attend to all possible binding sites; it improved blind-docking performance over EquiBind.
- **DiffDock** is a diffusion generative model over the non-Euclidean manifold of ligand poses; it outperformed previous methods on PDBbind and provides confidence estimates.
- **DynamicBind** uses equivariant geometric diffusion to generate a smooth energy landscape and unifies protein-conformation generation with ligand-pose prediction.

These methods are fast, but their physical interpretability and generalization remain under evaluation; see [[machine-learning-scoring-functions]].

## Practical algorithm choices

- For a rigid receptor and a drug-like ligand, incremental construction or a GA/MC search with a fast scoring function is standard.
- For a flexible pocket, ensemble docking is the most common practical compromise; soft docking is an alternative when only small side-chain rearrangements are expected.
- For blind docking, the search box covers the whole protein and the computational cost rises sharply; pocket prediction is usually preferable (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- AutoDock Vina uses a stochastic global search whose effort is controlled by the `exhaustiveness` parameter; the default is 8, and 32 is recommended for difficult ligands such as imatinib (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

## Related pages
- [[protein-ligand-docking]]
- [[scoring-functions]]
- [[autodock-vina]]
- [[binding-site-detection]]
- [[molecular-recognition-models]]
- [[molecular-dynamics-simulation]]
- [[machine-learning-scoring-functions]]
