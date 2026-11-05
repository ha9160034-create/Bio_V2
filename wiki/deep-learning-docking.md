# Deep Learning Docking

**Summary**: An overview of deep-learning generative, regression, and co-folding docking architectures, detailing their underlying mechanics, PoseBusters physical validity limits, and benchmark illusions.
**Sources**: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf; Integrating Machine Learning into Free Energy Perturbation Workflows.pdf; Protein–Ligand Docking in the Machine-Learning Era.md
**Last updated**: 2026-11-01

---

## Evolution of deep learning in molecular docking

Deep learning (DL) has emerged as an alternative to classical physics-based docking engines (such as Glide SP and AutoDock Vina), attempting to predict protein–ligand complex geometries and binding affinities directly from molecular structures without exhaustive physics-based conformational search (source: Protein–Ligand Docking in the Machine-Learning Era.md; Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf). Early machine-learning docking workflows focused on post-docking rescoring using 3D convolutional neural networks (e.g., Pafnucy, Gnina) and graph neural networks to extract non-linear spatial descriptors from docked poses (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Subsequent developments have introduced three distinct end-to-end algorithmic paradigms (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf):

### 1. Regression-based direct coordinate predictors
Early models (e.g., EquiBind, TANKBind, FABind) treat docking as a direct geometric regression task:
- EquiBind uses equivariant neural networks ($E(3)$-equivariance) to predict ligand center-of-mass translation, rotation, and internal torsion adjustments directly onto a rigid receptor pocket.
- TANKBind partitions protein pockets into sub-cavities and optimizes interatomic distance matrices using trigonometry-aware graph neural networks (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf).

### 2. Generative diffusion models
Generative architectures (such as DiffDock, DiffDock-L, DynamicBind, and DiffBindFR) frame docking as a reverse diffusion process over the ligand's SE(3) transformation manifold (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf):
- DiffDock performs diffusion over translational, rotational, and torsional degrees of freedom, generating candidate poses that are ranked by a separate confidence model.
- DynamicBind and NeuralPLexer introduce simultaneous pocket side-chain adjustments to model induced-fit transitions upon ligand insertion (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

### 3. All-atom biomolecular co-folding foundation models
Modern foundation models (such as AlphaFold 3, Chai-1, Boltz-1, and RoseTTAFold-All-Atom [RFAA]) unify structure prediction across proteins, nucleic acids, and small-molecule ligands:
- Model full atomistic interactions directly from primary amino acid sequences, multiple sequence alignments (MSAs), and ligand SMILES strings via unified pair-former and diffusion trunk architectures (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Integrating Machine Learning into Free Energy Perturbation Workflows.pdf).

## The PoseBusters benchmark and physical validity crisis

Large-scale evaluations by Morehead et al. (2026) and Li et al. (2025) demonstrate that reporting docking performance purely through root-mean-square deviation ($\text{RMSD} \le 2.0\text{ \AA}$) creates an exaggerated illusion of success (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf):
- **PoseBusters validity (PB-Valid)**: Measures whether predicted poses obey fundamental physical and chemical laws, testing for steric clashes, planar aromatic rings, non-planar aliphatic rings, valid covalent bond lengths, tetrahedral stereocenters, and torsional strain (source: Assessing the potential of deep learning for protein–ligand docking.pdf).
- **Physical validity collapse**:
  - KarmaDock achieved **0.00% PB-valid rates** on Astex Diverse and DockGen, and **0.23%** on PoseBusters, resulting in a **0.00% combined success rate** ($\text{RMSD} \le 2\text{ \AA}$ + PB-valid) across all benchmark datasets (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).
  - QuickBind yielded PB-valid rates of only 1.18% (Astex) and 0.00% (DockGen) (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).
  - Generative diffusion models (such as DiffDock) frequently place ligand atoms directly overlapping with receptor side chains (severe steric clashes) or generate unphysical ring conformations and strained bond lengths, requiring downstream force-field energy minimization (e.g., via OpenMM or Vina local minimization) to achieve PB-validity (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf).
- **Physics-based contrast**: In sharp contrast, conventional docking engines (AutoDock Vina, Glide SP) explore internal coordinate trees within empirical force fields, maintaining **100% PB-valid compliance** and often outperforming deep learning models on combined success metrics (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf).

## Generalization bottlenecks: memorization and the cross-docking gap

Critical benchmarks have exposed fundamental generalization limits in deep learning docking models (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf):

### 1. Dataset memorization and unseen pockets (DockGen)
When tested on non-redundant, unseen protein receptor folds (DockGen dataset), DL docking performance drops steeply:
- SurfDock's combined success rate falls from 61.18% (Astex) to 33.33% (DockGen);
- Interformer-Energy drops from 68.24% to 34.39%;
- DiffBindFR collapses from 65.88% to 18.52% (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).
These steep drops reveal that DL models memorize ligand-pocket co-occurrences in the PDBbind training set rather than learning generalizable biophysical interaction potentials (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).

### 2. The redocking versus cross-docking and apo-docking collapse
Almost all DL docking architectures were trained and benchmarked on **holo crystal structures**, where the protein binding site is already pre-formed to fit the co-crystallized ligand ("self-docking" or "redocking") (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf):
- When evaluated under realistic **cross-docking** (docking ligands into alternate holo receptor structures) or **apo-docking** (docking into unliganded, apo crystal structures), DL success rates collapse (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).
- DL models lack explicit force-field repulsion mechanics to guide the clearance of unrelaxed side-chain rotamers, leading to unresolvable steric clashes (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).

### 3. Evolutionary MSA dependency in co-folding models
Morehead et al. (2026) demonstrated that co-folding foundation models (AlphaFold 3, Chai-1, Boltz-1) rely heavily on deep Multiple Sequence Alignments (MSAs). When evaluated with single-sequence inputs (without evolutionary homologues), pose accuracy plummets, proving that the models infer binding pocket architecture from evolutionary co-conservation rather than de novo physical chemical forces (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

### 4. Virtual screening failure
In structure-based virtual screening campaigns designed to enrich active drugs over property-matched decoys:
- DL docking scores show weak or non-existent correlation with experimental binding affinities;
- Models are easily fooled by decoy molecules with high atomic surface contact, often performing **no better than random guessing** and failing to match conventional physics-based scoring functions (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).

## Related pages
- [[protein-ligand-docking]]
- [[machine-learning-scoring-functions]]
- [[benchmarking-datasets]]
- [[docking-failure-modes]]
- [[autodock-vina]]
- [[virtual-screening]]
- [[alchemical-free-energy-methods]]
- [[morehead-2026-evaluating-dl-docking]]
- [[li-2025-decoding-dl-docking-limits]]
