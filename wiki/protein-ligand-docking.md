# Protein–Ligand Docking

**Summary**: Molecular docking predicts the non-covalent binding pose and approximate affinity of a small molecule in a protein binding site by combining conformational search with a scoring function.
**Sources**: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Protein–Ligand Docking in the Machine-Learning Era.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Basic docking — Autodock Vina 1.2.0 documentation.md; Molecular dynamics simulations and drug discovery.md; Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf; The Light and Dark Sides of Virtual Screening.pdf; Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf; Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf
**Last updated**: 2026-11-01

---

## Definition

Protein–ligand docking is a computational procedure that predicts the non-covalent complex formed between a macromolecular receptor and a small molecule. It is described as solving a three-dimensional jigsaw puzzle: the algorithm searches for the association state that best matches the two molecules and then estimates how tightly they bind (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Docking is a central tool in computer-aided drug design (CADD). It is used for binding-mode prediction, structure–activity and mutagenesis studies, virtual screening, and lead optimization (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## The two essential components

Every docking calculation reduces to **sampling** and **ranking** (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md):

1. **Conformational sampling** generates putative ligand positions, orientations, and torsion angles inside the binding site. A ligand has three translational and three rotational degrees of freedom, plus one degree of freedom for each rotatable dihedral; proteins add far more flexibility (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
2. **Scoring** ranks the generated poses and estimates binding affinity. The scoring function is the component that determines whether a docking program can distinguish a native pose from a decoy (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

Most docking programs keep the receptor rigid and treat only the ligand as flexible, because a protein contains thousands of atoms and full receptor sampling is computationally prohibitive (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). Receptor flexibility is handled indirectly by ensemble docking, soft docking, or by restricting flexibility to selected side chains; see [[molecular-docking-algorithms]] and [[virtual-screening]].

## What docking predicts

A docking run produces a ranked list of ligand **poses** (typically written as PDBQT files containing atomic coordinates and partial charges) together with an estimated binding affinity in kcal/mol (source: Basic docking — Autodock Vina 1.2.0 documentation.md). In standard implementations like AutoDock Vina, the docking engine searches over ligand translational, rotational, and active dihedral torsion degrees of freedom within a pre-defined 3D search box, using iterated local search or genetic algorithms to optimize the score (source: Basic docking — Autodock Vina 1.2.0 documentation.md). The resulting score is not a true free energy: it is a fast approximation built from simplified interaction terms and usually omits or crudely treats solvent and entropy (source: Basic docking — Autodock Vina 1.2.0 documentation.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

The three classical applications of a scoring function in docking are (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md):

- determining the binding mode and binding site of a ligand;
- predicting absolute binding affinity, chiefly for lead optimization;
- ranking a compound database in structure-based [[virtual-screening]].

## Accuracy: sampling versus scoring and self-docking vs blind challenges

The apparent accuracy of molecular docking depends heavily on benchmark design:
- **Self-docking / Redocking benchmarks**: In retrospective self-docking where a co-crystallized ligand is redocked into its pre-formed holo receptor cavity, physics-based programs routinely achieve 70% to 85% success rates ($\text{RMSD} \le 2.0\text{ \AA}$) across standard test sets (such as Astex Diverse) (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- **Prospective blind challenges (CELPP)**: In the Continuous Evaluation of Ligand Pose Prediction (CELPP) community exercise, the Zou group evaluated blind docking on weekly newly deposited, unpublished PDB crystal structures (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). CELPP revealed a major performance dichotomy: **conformational sampling succeeded in 80.4% of cases, whereas scoring succeeded in only 46.3%** of targets (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). In other words, while algorithms reliably generate near-native poses within their search ensembles, scoring functions frequently misrank decoys above the true binding mode.

The authors attribute scoring errors to solvent displacement effects, entropic penalties, ambiguous protonation/tautomeric states, and neglected electronic polarization (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). They identify **protein flexibility as the biggest remaining hurdle**, because sampling all relevant receptor conformations is computationally immense (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

A separate benchmark of scoring functions on 100 diverse complexes illustrates how uneven performance can be: ITScore/SE ranked 91% of native poses correctly (RMSD ≤ 2.0 Å), whereas SYBYL/D-Score reached only 26%; for affinity correlation, ITScore/SE and ITScore reached R = 0.65, while AutoDock reached R = 0.05 on that particular set (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md). These values are test-set dependent and should not be read as universal rankings; see [[scoring-functions]] and [[benchmarking-datasets]].

## Docking in the drug-discovery workflow

Docking sits between experimental structural biology and wet-lab screening. It requires a three-dimensional receptor structure, which may come from X-ray crystallography, cryo-EM, NMR, homology modelling, or deep-learning structure prediction tools such as AlphaFold2/RoseTTAFold and modern all-atom co-folding models (AlphaFold 3, Chai-1, Boltz-1) (source: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

A typical structure-based workflow is: obtain and prepare the receptor; identify or predict the binding site; prepare ligands with correct protonation and geometry; dock with a flexible-ligand/rigid-receptor protocol, often across an ensemble of receptor conformations; then rescore or filter the poses before experimental testing (source: Protein–Ligand Docking in the Machine-Learning Era.md). The full workflow is described in [[virtual-screening]] and [[structure-based-drug-design]].

## Relationship to molecular dynamics and free-energy methods

Docking is optimized for speed: it evaluates a single protein–ligand complex with a simplified score. Molecular dynamics (MD) and free-energy calculations sample the actual trajectory and solvent, and are therefore more accurate but orders of magnitude more expensive (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Molecular dynamics simulations and drug discovery.md; Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf). Bridging docking with MD simulations provides vital validation of docked poses, resolves ambiguous electron density, and tracks the kinetic pathway of ligand recognition and entry (source: Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf).

The practical hierarchy is:

- **Docking + scoring function**: microseconds to milliseconds per ligand; used for ultra-large libraries; accuracy limited by the score (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- **MM-PBSA / MM-GBSA**: post-processing on MD snapshots; intermediate accuracy and cost; see [[mm-pbsa-and-mm-gbsa]].
- **FEP / TI / alchemical perturbation**: rigorous free-energy calculations; highest cost, used for congeneric series in lead optimization; see [[alchemical-free-energy-methods]].
- **PMF / umbrella sampling / metadynamics**: physically meaningful pathways that also yield kinetics; see [[enhanced-sampling]].

These methods are detailed in [[binding-free-energy]], [[molecular-dynamics-simulation]], and [[enhanced-sampling]].

## Physics-based versus deep-learning docking

Recent evaluations benchmark traditional physics-based engines (Glide SP, AutoDock Vina, DOCK 3.7) against emerging deep learning docking and co-folding algorithms (DiffDock, Chai-1, Boltz-1, AlphaFold 3) (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf):
- **PoseBusters physical validity**: Generative DL models frequently generate severe steric clashes, impossible bond lengths, and distorted rings (source: Assessing the potential of deep learning for protein–ligand docking.pdf). In large-scale benchmarks, KarmaDock achieved **0.00% PB-valid compliance** across Astex Diverse and DockGen datasets, yielding 0.00% combined success (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf). DiffDock places ligand atoms directly overlapping with receptor side chains in a high fraction of poses, requiring downstream force-field energy minimization to achieve validity (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf). In contrast, physics-based engines (AutoDock Vina, Glide SP) maintain near 100% geometric validity.
- **Cross-docking and apo-docking collapse**: DL models trained on holo crystal structures collapse when presented with apo receptors or alternate holo conformations, as they lack explicit physical force fields to clear unrelaxed side-chain rotamers (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf). Detailed analysis is provided in [[deep-learning-docking]].

## Key failure modes and limitations

Systematic benchmarks delineate the recurring failure modes of molecular docking (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf; The Light and Dark Sides of Virtual Screening.pdf):
- **Protein flexibility & loop occlusion**: Rigid-receptor docking fails when binding requires side-chain clearance or loop displacement (e.g., MAPK2 80s loop) (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Empirical molecular weight bias**: Vina's empirical scoring function strongly correlates with molecular weight, systematically penalizing low-MW fragments and rewarding large compounds (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Ligand strain & irregular torsions**: Docking programs that omit calibrated internal strain energy penalties generate distorted dihedrals (311 and 281 irregular torsion types in Vina and DOCK 3.7) (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Metalloenzyme coordination errors**: Standard docking protocols zero out or misassign metal charges, causing complete active enrichment failure on targets like COMT unless explicit partial charges or specialized metal potentials are supplied (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Chemical artifacts**: False positives from Pan-Assay Interference Compounds (PAINS) and colloidal aggregators (source: The Light and Dark Sides of Virtual Screening.pdf). See [[docking-failure-modes]].

## Related pages
- [[molecular-docking-algorithms]]
- [[scoring-functions]]
- [[deep-learning-docking]]
- [[docking-failure-modes]]
- [[mm-pbsa-and-mm-gbsa]]
- [[alchemical-free-energy-methods]]
- [[machine-learning-scoring-functions]]
- [[autodock-vina]]
- [[virtual-screening]]
- [[structure-based-drug-design]]
- [[binding-site-detection]]
- [[protein-ligand-interactions]]
- [[binding-free-energy]]
- [[molecular-dynamics-simulation]]
- [[benchmark-drug-targets]]
- [[salmaso-2018-bridging-docking-md]]
- [[ma-2024-docking-strategies]]
