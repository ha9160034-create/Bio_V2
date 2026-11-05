# Benchmarking Datasets

**Summary**: Scoring functions and docking protocols are evaluated on curated protein–ligand datasets that pair 3D structures with experimental affinities or with experimentally verified active/inactive labels; PDBbind/CASF, DUD-E, MUV, LIT-PCBA, CSAR, D3R, and CACHE cover different evaluation tasks and carry different biases.
**Sources**: Protein–Ligand Docking in the Machine-Learning Era.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf; Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf; The Light and Dark Sides of Virtual Screening.pdf; Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf; Assessment of molecular dynamics time series.pdf; The Open Force Field Initiative Open Software and Open Science.pdf; Training a force field for proteins and small molecules from.pdf
**Last updated**: 2026-11-02

---

## The three evaluation tasks

A scoring function is judged on three related but distinct tasks (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Protein–Ligand Docking in the Machine-Learning Era.md):

1. **Binding pose identification** — can the native pose be ranked first among decoys?
2. **Binding affinity prediction** — does the score correlate with experimental affinity?
3. **Virtual screening** — can true actives be enriched above inactives?

A function that excels at one task may fail at another. RF-Score, for example, outperformed classical scoring functions on affinity prediction but performed much worse on docking and screening (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Structure–affinity datasets

### PDBbind

PDBbind is the largest curated collection of 3D protein–ligand structures with experimentally measured binding data (K_d, K_i, or IC50), collected manually from the original references and updated annually (source: Protein–Ligand Docking in the Machine-Learning Era.md).

- **General set, version 2020:** 19,443 complexes with affinities from **1.2 pM to 10 mM** (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Refined set:** a high-quality subset selected by structure and affinity criteria (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Core set:** the basis of the Comparative Assessment of Scoring Functions (CASF) (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### CASF-2016

CASF-2016 is derived from the PDBbind refined set and contains **57 targets with 5 crystal complexes each** — 285 complexes total, with ligands spanning high to low affinity for each target (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Scoring functions are evaluated on four metrics:

- **scoring** — Pearson correlation between predicted and experimental affinity (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- **ranking** — average Spearman rank correlation across the 57 targets (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- **docking** — fraction of predicted poses within 2 Å RMSD of the crystal pose (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- **screening** — success in ranking the highest-affinity binder first, plus average enrichment factor at 1%, 5%, and 10% (source: Protein–Ligand Docking in the Machine-Learning Era.md).

The screening assessment is limited by the small database (285 compounds) and the lack of experimentally verified inactives (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### CSAR

The Community Structure–Activity Resource provides exercises and high-quality affinity data for validating scoring functions and SBDD methods (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### D3R Grand Challenge

The Drug Design Data Resource ran blind, open competitions from 2015 to 2019 for pose prediction and affinity ranking. Grand Challenge 4 covered macrocyclic BACE1 inhibitors (pose prediction in stages 1a and 1b, affinity ranking in stage 2) and Cathepsin S inhibitors (affinity ranking). Mean and median RMSD assess poses; Spearman and Kendall correlations assess affinity ranking (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### CACHE

The Critical Assessment of Computational Hit-finding Experiments launched in early 2022 as a prospective benchmark for virtual screening in real screening campaigns. It organizes multiple rounds against novel targets, with experimental testing of predictions (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### PoseBusters Benchmark
Developed by Buttenschoen et al., PoseBusters evaluates the chemical and physical plausibility of docked protein–ligand complexes (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf):
- Consists of curated test sets ($n = 130$ or $308$ complexes) filtered to exclude crystallography artifacts;
- Imposes **PB-Valid** criteria: checks for steric clashes, planar aromatic rings, non-planar aliphatic rings, valid covalent bond lengths, tetrahedral stereocenters, and torsional strain;
- Establishes that generative deep learning models often achieve low RMSD while failing PB-Valid checks (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

### Astex Diverse Set
A high-resolution, structurally diverse benchmark of 85 structurally distinct protein–ligand complexes widely adopted as a standard test set for crystallographic redocking and pose reproduction (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf).

### DockGen and DockGen-E
Curated by Corso et al. to test the out-of-distribution generalization of machine learning models:
- Enforces structural partitioning where test receptor folds share no sequence or structural homology with the PDBbind training set;
- Exposes the collapse of deep-learning models when confronted with novel, unmemorized binding pockets (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).

### PLBenchmarks and arsenic
Developed by Hahn et al. (2022), **PLBenchmarks** provides an open-access, versioned, standardized benchmark repository for evaluating relative (RBFE) and absolute (ABFE) alchemical binding free energy calculations (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf):
- Filters experimental bioactivity strictly to direct biophysical assays (ITC, SPR) and high-quality binding assays, rejecting uncalibrated enzymatic IC50 values (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf);
- Enforces an experimental affinity dynamic range spanning $> 2\text{--}3$ orders of magnitude ($\Delta\Delta G \ge 2.7\text{--}4.1\text{ kcal/mol}$), keeping in mind the experimental reproducibility noise floor of $\sim 0.3\text{--}0.5\text{ kcal/mol}$ (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf);
- Requires high-resolution crystal structures ($\le 2.0\text{--}2.5\text{ \AA}$, {\text{free}} < 0.25$, RSCC density validation) and closed thermodynamic perturbation cycles (cycle closure hysteresis $\sum_{\text{cycle}} \Delta\Delta G = 0$) (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf);
- Pairs with **arsenic**, an open-source evaluation toolkit calculating RMSE, MAE, Kendall $\tau$, and bootstrapped 95% confidence intervals while avoiding misleading ^2$ inflation (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). See [[simulation-and-benchmark-best-practices]] and [[hahn-2022-binding-affinity-benchmarks]].

### The 2,752-Complex MD Time-Series Benchmark
Curated by Poziemski et al. (2026), this dataset encompasses explicit-solvent molecular dynamics simulations across **2,752 protein–ligand complexes** from PDBbind to evaluate dynamical machine learning scoring functions (source: Assessment of molecular dynamics time series.pdf):
- Extracts dynamic time-series features (	s_descriptors) using 	sfresh across 75 biophysical metrics, filtering out simulation noise via p-value trimming, correlation pruning (PCC $\ge 0.8$), and mutual information maximization (source: Assessment of molecular dynamics time series.pdf);
- Implements rigorous **Target Splits** (UniProt accession IDs) and **Scaffold Splits** (Bemis–Murcko chemotype clusters via DeepChem) to benchmark out-of-distribution generalizability on unseen biological families and chemical series (source: Assessment of molecular dynamics time series.pdf). See [[poziemski-2026-md-time-series-descriptors]].

### FreeSolv and SPICE Thermophysical Benchmarks
Used extensively by OpenFF and Garnet for force field validation:
- **FreeSolv**: A curated database of 642 neutral small molecules with experimentally measured aqueous hydration free energies, used to benchmark force-field electrostatics and solvation parameters (source: The Open Force Field Initiative Open Software and Open Science.pdf);
- **SPICE (Secondary Pseudopotential Interatomic Coordinate Energy)**: A large-scale quantum chemical dataset comprising hundreds of thousands of DFT single-point energies, forces, and conformational profiles across drug-like small molecules, dipeptides, and molecular dimers, used to train continuous neural network force fields like Garnet and benchmark valence parameters in OpenFF Sage (source: Training a force field for proteins and small molecules from.pdf). See [[open-force-field-and-ligand-parameterization]].

### CASP15 PLI
The Critical Assessment of Structure Prediction (CASP15) Protein–Ligand Interface challenge tests blind prediction of complexes where target structures and bound ligands are withheld until after computational submission (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

## Active/inactive datasets

### DUD and DUD-E

The Database of Useful Decoys and its Enhanced version are widely used for VS benchmarking. DUD-E contains **102 targets**, **22,886 active compounds**, and, for each active, **50 computer-generated decoys** with similar physicochemical properties but dissimilar 2D topology (source: Protein–Ligand Docking in the Machine-Learning Era.md; Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

The major drawback is that many decoys are presumed inactive without experimental verification, so false negatives may exist (source: The Light and Dark Sides of Virtual Screening.pdf). Decoys also differ systematically from actives in topology and charge distribution, which can artificially inflate apparent enrichment factors (source: Protein–Ligand Docking in the Machine-Learning Era.md; The Light and Dark Sides of Virtual Screening.pdf). In a comprehensive evaluation across all 102 DUD-E targets, Xu et al. revealed that target pocket rigidity, loop occlusions, and molecular weight biases account for substantial failure rates across diverse target classes (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf). See [[docking-failure-modes]] and [[benchmark-drug-targets]].

### MUV

Maximum Unbiased Validation was built from PubChem bioactivity data for **17 targets**, each with **30 actives** and **15,000 inactives**, and was designed to avoid analogue bias and artificial enrichment. Unlike DUD/DUD-E, MUV provides experimentally verified inactives, mostly from cell-based assays. However, many actives are not validated against their putative targets, so MUV is more appropriate for ligand-based VS than structure-based VS (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### LIT-PCBA

LIT-PCBA was curated in 2020 from PubChem dose–response assays. It contains **15 targets**, with all actives and inactives taken from experimental data under homogeneous conditions. Potential false positives were carefully removed by requiring dose–response curves with Hill slopes between 0.5 and 2.0 (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Its main limitation is that **8 of the 15 targets are cell-based phenotypic assays**, which introduces membrane permeability, cytotoxicity, and other whole-cell effects into a structure-based benchmark (source: Protein–Ligand Docking in the Machine-Learning Era.md).

LIT-PCBA results discussed in the wiki (source: Protein–Ligand Docking in the Machine-Learning Era.md):

- IFP post-docking filtering: average **EF1% = 7.46**, the best early enrichment (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- GRIM interaction graph matching: similar performance (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- ΔLin_F9 XGB: **13 of 15 targets with EF1% > 2**, average **EF1% = 5.55**, best among affinity-predicting methods (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- template-based methods (FINDSITEcomb2.0, FragSite) and GNINA achieved lower enrichment (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Comparisons are complicated by different docking protocols — conformer generation, box definition, and number of templates — which can change Vina-based EF1% between studies (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Large bioactivity databases

These contain many compounds and affinities but few or no annotated 3D structures (source: Protein–Ligand Docking in the Machine-Learning Era.md):

- **BindingDB** — as of 4 May 2022, 2,513,948 binding data for 8,839 protein targets and 1,077,922 small molecules (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **ChEMBL** — version 30 contains 19,286,751 activities for 14,885 targets and 2,157,379 compounds (source: Protein–Ligand Docking in the Machine-Learning Era.md).

They are used for ligand-based and sequence-based affinity prediction and can supplement structure-based scoring-function development.

## Evaluation metrics

### Retrospective redocking versus blind prospective scoring (CELPP)

Benchmark performance diverges dramatically between retrospective crystallographic redocking and blind prospective scoring (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf):
- **Self-docking / Redocking**: Docking a ligand back into its native, pre-formed holo crystal structure achieves high nominal success rates (typically 70% to 85% of complexes with RMSD $\le 2.0\text{ \AA}$) across standard test sets (such as Astex Diverse or PDBbind core sets) (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; The Light and Dark Sides of Virtual Screening.pdf).
- **Prospective blind challenges (CELPP)**: The Continuous Evaluation of Ligand Protein Predictions (CELPP) community challenge assessed docking performance on weekly newly deposited PDB structures before publication (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf). CELPP demonstrated that while conformational search algorithms successfully sampled near-native poses ($\text{RMSD} \le 2.0\text{ \AA}$) in **80.4%** of cases, existing scoring functions identified the correct pose as top-ranked in only **46.3%** of targets (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf). This underscores that scoring and ranking—rather than conformational sampling—constitute the primary bottleneck in prospective structure-based design (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf; The Light and Dark Sides of Virtual Screening.pdf).

### Pose prediction

The root-mean-square deviation between the top docking pose and the experimental pose. A prediction is successful if **RMSD ≤ 2 Å** (source: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

For symmetric molecules, the **minimum symmetry-corrected RMSD** over all equivalent atom mappings must be used; otherwise a correct pose can be scored as wrong (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### Affinity prediction

- **Pearson correlation coefficient** between predicted scores and experimental affinities — assumes a linear relationship.
- **Spearman rank correlation** — ranks both sets first and correlates the ranks; appropriate when the relationship is monotonic but not linear.

Some scoring functions produce scores that cannot be compared directly with experimental binding data, in which case correlation, not absolute error, is the appropriate metric (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### Virtual screening

- **Enrichment factor (EF)** — the fraction of true binders found in the top α percentile divided by the fraction expected at random:

$$EF_\alpha = \frac{N_{actives,\ \alpha}/N_{actives}}{N_{total,\ \alpha}/N_{total}}$$

- **ROC curve and AUC** — true-positive rate versus false-positive rate; AUC = 0.5 is random, 1.0 is perfect. AUC is most appropriate when actives and inactives are comparable in number (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Dataset construction criteria

A good training or test set should satisfy (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md):

- high-quality structures with no atomic clashes, preferably high-resolution crystal structures;
- a wide range of proteins and binding affinities;
- drug-like, non-covalent ligands;
- experimentally verified actives and inactives;
- no overlap between training and test sets;
- a composition that reflects the intended screening library.

## Bias and pitfalls

- **Decoy bias** — DUD-E decoys are presumed inactive and topologically dissimilar; enrichment can be inflated (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **False negatives** — unverified inactives may include true binders.
- **Analogue bias** — actives from the same scaffold cluster; MUV was designed to avoid this (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Assay heterogeneity** — IC50 values from different assays, cell lines, and conditions are not equivalent to K_d.
- **Cell-based confounds** — LIT-PCBA includes phenotypic assays where permeability and toxicity affect the label (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Small benchmark size** — CASF's 285 compounds limit screening assessment (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Protocol dependence** — enrichment depends on docking box, conformer generation, and receptor structure (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Train/test leakage** — homologous targets or overlapping ligands can make ML scores look better than they are on novel targets.
- **Retrospective versus prospective** — most results are retrospective; CACHE and D3R provide prospective tests (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## How to report a benchmark

1. Name the dataset and version (source: Protein–Ligand Docking in the Machine-Learning Era.md).
2. State the train/test split and any overlap (source: Protein–Ligand Docking in the Machine-Learning Era.md).
3. Describe receptor preparation, ligand preparation, box definition, and docking program (source: Protein–Ligand Docking in the Machine-Learning Era.md).
4. Report pose success (RMSD ≤ 2 Å, symmetry-corrected), affinity correlation (Pearson and Spearman), and enrichment (EF1%, EF5%, AUC) (source: Protein–Ligand Docking in the Machine-Learning Era.md).
5. Include confidence intervals or standard deviations across targets (source: Protein–Ligand Docking in the Machine-Learning Era.md).
6. Report failures, not only successes (source: Protein–Ligand Docking in the Machine-Learning Era.md).
7. Validate prospectively where possible (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Related pages
- [[scoring-functions]]
- [[machine-learning-scoring-functions]]
- [[deep-learning-docking]]
- [[docking-failure-modes]]
- [[virtual-screening]]
- [[protein-ligand-docking]]
- [[experimental-binding-assays]]
- [[benchmark-drug-targets]]
- [[yang-2022-ml-docking]]
- [[huang-2010-scoring-functions]]
- [[morehead-2026-evaluating-dl-docking]]
- [[li-2025-decoding-dl-docking-limits]]
