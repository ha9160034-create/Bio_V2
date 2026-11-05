# Machine-Learning Scoring Functions

**Summary**: Machine-learning scoring functions replace or correct classical linear scoring terms with learned models ranging from random forests to graph neural networks, achieving strong affinity prediction but uneven docking and screening performance.
**Sources**: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf; Integrating Machine Learning into Free Energy Perturbation Workflows.pdf; Assessment of molecular dynamics time series.pdf; Training a force field for proteins and small molecules from.pdf
**Last updated**: 2026-11-02

---

## From linear scores to learned scores

Classical scoring functions — physics-based, knowledge-based, and empirical — usually adopt a linear functional form: a weighted combination of force-field terms, contact counts, or empirical interaction descriptors (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Protein–Ligand Docking in the Machine-Learning Era.md). While linear functions offer computational efficiency and clear physical interpretability, they struggle to model non-linear cooperativity, multi-body desolvation effects, and complex electronic polarization (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md). Machine-learning (ML) scoring functions overcome these rigid functional constraints using non-linear models ranging from support vector machines and random forests to deep convolutional networks and graph neural networks (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Protein–Ligand Docking in the Machine-Learning Era.md).

The key promise is that ML can learn corrections to an existing physical score rather than replacing it entirely. The key risk is that a model trained on crystal structures may learn dataset artefacts rather than transferable binding physics (source: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## Early descriptor-based models

- **RF-Score** (Ballester and Mitchell, 2010) was the first ML scoring function to outperform classical scoring functions on the scoring task. It used random forests with protein–ligand atom-type pair counts within distance cutoffs, trained on PDBbind v2007 (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **SFCscore^RF** (2013) used 63 empirical descriptors covering ligand properties, specific interactions such as hydrogen bonds and aromatic contacts, and surface characteristics. It outperformed classical functions on CASF 2013 and CSAR-NRC HiQ (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **NNScore 1.0** used a single hidden layer with five neurons and 194 features to classify active versus inactive compounds. **NNScore 2.0** used more interaction terms and predicted pK_d rather than a binary class (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Both RF-Score and SFCscore^RF performed much worse on docking and screening than on affinity prediction (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Δ-machine learning: correcting Vina

The Zhang group introduced a Δ-machine-learning strategy in which an ML model parameterizes corrections to the Vina score. This retains Vina's docking power while improving affinity prediction (source: Protein–Ligand Docking in the Machine-Learning Era.md).

- **ΔVina RF20** (2017) used 10 pharmacophore-based solvent-accessible surface area (SASA) features and 10 of Vina's 58 empirical terms. It achieved the best performance among classical scoring functions across scoring, ranking, docking, and screening on CASF 2007 and 2013.
- **ΔVina XGB** (2019) replaced random forest with extreme gradient boosting and added explicit water molecules, ligand-conformation stability, and metal terms. Its feature set comprised Vina's 58 features, 30 SASA features, 3 bridge-water features, 2 ligand-stability features, and 1 metal-count term. Training included explicitly solvated and dry structures plus docking decoys; it consistently improved all four CASF-2016 metrics.
- **ΔLin_F9 XGB** (2022) used Gaussian terms to characterize interactions at different distance ranges and enlarged the training set with weak binders and docking poses. It achieved superior scoring, ranking, and screening on CASF-2016.

## Interaction-feature models

- **ET-score** (2021) used distance-weighted interatomic contacts between protein and ligand atom-type pairs and reached Pearson R = 0.827 on CASF-2016.
- **ECIF-GBT** (2021) used extended connectivity interaction features — atom-type pair counts that account for each atom's connectivity — and reached Pearson R = 0.857 on CASF-2016.

Both were trained solely on crystal structures, and their docking and screening performance remains an issue (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Deep learning: grids and graphs

### Convolutional models

- **AtomNet** (2017) was the first CNN-based scoring function to incorporate 3D structural information. It used vectorized 3D grids over the protein–ligand interface, with four 3D convolutional layers and two fully connected layers, and achieved a better AUC than Vina on DUD-E.
- **Pafnucy** (2017) and **Kdeep** (2018) used atom-property or atom-type 3D grids. Kdeep maintained the 3D structure of the input grid and embedded more physical information than AtomNet; its accuracy was competitive with other scoring methods but significantly below FEP.
- A limitation of grid CNNs is dependence on the coordinate frame: different orientations of the same complex can produce different representations.

### Rotation-invariant and graph models

- **OnionNet** (2019) used rotation-free element-pair-specific contacts in different distance shells; OnionNet-2.0 followed in 2021 and performed well on CASF-2016.
- **PotentialNet** (2018) used a graph convolutional network to learn intramolecular and intermolecular interactions through covalent-only, dual non-covalent/covalent, and ligand-graph-gathering stages; it reached Pearson R = 0.822 on CASF-2007.
- **graphDelta** (2020) and **SIGN** (2021) are later graph models; however, these GNNs did not outperform some traditional ML scoring functions on CASF benchmarks (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **RTMScore** used a residue–atom distance likelihood potential with a graph transformer and performed strongly on CASF-2016 docking and screening (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## The generalization problem

A central critique is that ML scoring functions can appear excellent on affinity benchmarks while failing at the tasks that matter for virtual screening. Gabel et al. reproduced the claimed scoring power of RF-Score and SVM-based functions but found their docking power (discriminating near-native from decoy poses) and screening power (discriminating actives from decoys) to be poor (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

The problem has several sources:

- **Dataset bias**: decoys may be physically unrealistic or easy to distinguish; DUD-E decoys are presumed inactive without experimental verification (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Training-set leakage**: closely related proteins or ligands can appear in both training and test sets.
- **Black-box behaviour**: deep models are difficult to interpret and can be over-trained (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Novel targets**: there is no guarantee that an ML function trained on PDBbind will outperform a classical function on a target unlike anything in its training data (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## From scoring functions to end-to-end deep learning docking and active learning

Beyond post-docking rescoring, deep learning has expanded in two major directions:
1. **Active learning surrogate scorers**: Machine-learning scoring functions (GNNs, random forests, Bayesian regression) act as fast surrogate filters coupled with rigorous alchemical Free Energy Perturbation (FEP) engines (source: Integrating Machine Learning into Free Energy Perturbation Workflows.pdf). By training on a subset of computed alchemical free energies and iteratively prioritizing compounds with high acquisition values, ML-FEP workflows achieve **20-fold speedups**, enabling alchemical accuracy across multi-million-compound libraries (source: Integrating Machine Learning into Free Energy Perturbation Workflows.pdf).
2. **Direct generative and co-folding docking**: Architectures such as DiffDock, KarmaDock, SurfDock, and AlphaFold 3 predict coordinates directly from sequence or molecular graphs (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf). However, rigorous multi-benchmark evaluations have uncovered critical physical and generalization limitations:
- **PoseBusters physical validity crisis**: Generative models frequently report high nominal pose accuracy (RMSD $\le 2.0$ Å) while generating severe steric clashes, distorted ring planes, and unphysical bond lengths (source: Assessing the potential of deep learning for protein–ligand docking.pdf). The PoseBusters benchmark (PB-Valid) evaluates whether poses satisfy physical chemistry rules (steric overlap ratio < 0.075, bond lengths within 4 standard deviations, correct aromatic planarity, and lack of internal torsional clash) (source: Assessing the potential of deep learning for protein–ligand docking.pdf).
- **Physical validity collapse in DL models**: KarmaDock achieved **0.00% PB-valid rates** on Astex Diverse and DockGen, resulting in a **0.00% combined success rate** ($\text{RMSD} \le 2\text{ \AA}$ + PB-valid) across benchmark datasets (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf). DiffDock frequently produces severe steric clashes between ligand atoms and binding pocket residues, requiring downstream force-field energy minimization to achieve physical validity (source: Assessing the potential of deep learning for protein–ligand docking.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf).
- **Out-of-distribution failure on unseen pockets (DockGen)**: Performance drops precipitously on novel folds (e.g., SurfDock combined success dropped from 61.18% to 33.33%), exposing heavy reliance on memorizing PDBbind complexes rather than learning biophysical forces (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).
- **Cross-docking and apo collapse**: DL models trained on holo crystal structures collapse when docking into unbound (apo) or alternate holo structures, as they lack physical force fields to clear unrelaxed side chains (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf). See [[deep-learning-docking]] and [[docking-failure-modes]].

## Case study: LIT-PCBA

The ΔLin_F9 XGB protocol was tested on the LIT-PCBA benchmark, which contains 15 targets with curated active and inactive compounds from PubChem dose–response assays (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Protocol: the co-crystal ligand defined the docking box; RDKit generated and protonated low-energy 3D conformers; flexible-ligand/rigid-receptor docking used Smina with the Lin_F9 scoring function; the top five poses were rescored with ΔLin_F9 XGB, and the best-rescored pose was used for assessment (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Results:

- **Interaction-fingerprint (IFP) filtering** achieved the best average early enrichment: EF1% = 7.46 across all 15 targets (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **GRIM** (interaction graph matching) performed similarly (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **ΔLin_F9 XGB** achieved average EF1% = 5.55 and the greatest number of targets with EF1% > 2 (13 of 15), representing the best method among those that predict binding affinity (source: Protein–Ligand Docking in the Machine-Learning Era.md).

The comparison is complicated by different docking protocols: IFP and GRIM depend on a reference crystal pose and do not predict absolute binding free energy, while the ML scoring approach does (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Molecular dynamics time-series descriptors (ts_descriptors)

A pervasive drawback of conventional ML scoring functions is their strict reliance on static crystal coordinates, which ignore thermal fluctuations, entropic interaction persistence, and induced-fit pocket breathing (source: Assessment of molecular dynamics time series.pdf).

Poziemski et al. (2026) systematically evaluated the impact of MD time-series descriptors (	s_descriptors) across **2,752 protein–ligand complexes** from PDBbind using Random Forest, XGBoost, and Multilayer Perceptrons (source: Assessment of molecular dynamics time series.pdf):
- **Simulation Noise vs. Curated Dynamics**: Feeding raw trajectory frames directly into ML models introduces severe simulation noise that degrades scoring accuracy as trajectory length increases (source: Assessment of molecular dynamics time series.pdf).
- **Five-Stage Feature Selection**: Extracting time-series features with 	sfresh across 75 biophysical metrics, followed by $-value trimming, correlation pruning (removing collinear pairs with Pearson correlation coefficient $\text{PCC} \ge 0.8$), and mutual information filtering, isolates compact vectors of at most 51 highly informative dynamic descriptors (source: Assessment of molecular dynamics time series.pdf).
- **Generalization Regimes**: While static crystal features perform well on random cross-validation splits, dynamic 	s_descriptors significantly improve generalization on challenging **Target Splits** (unseen protein families) and **Scaffold Splits** (unseen ligand chemotypes) by capturing physical interaction stability over time (source: Assessment of molecular dynamics time series.pdf). See [[poziemski-2026-md-time-series-descriptors]].

In parallel, machine learning has expanded to continuous force-field parameterization: **Garnet** uses graph neural networks and double exponential nonbonded potentials to learn all molecular mechanics parameters from scratch for both proteins and small molecules, bridging machine learning directly into rigorous physical MD simulation (source: Training a force field for proteins and small molecules from.pdf). See [[blanco-gonzalez-2024-garnet-force-field]].

## Target-specific models

Generic ML scoring functions aim to work on all targets, which is difficult given the particularities of each binding site. Target-specific ML functions trained on a single target or family can achieve state-of-the-art performance when sufficient data exist, but they are not applicable to novel targets with little experimental data (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Practical assessment

ML scoring functions are most useful as **rescorers** of poses generated by a physics-based docking program, especially when the training data resemble the target family. They should be validated on docking and screening metrics, not only on affinity correlation, and their predictions on novel chemotypes should be treated with caution (source: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## Related pages
- [[scoring-functions]]
- [[poziemski-2026-md-time-series-descriptors]]
- [[blanco-gonzalez-2024-garnet-force-field]]
- [[deep-learning-docking]]
- [[simulation-and-benchmark-best-practices]]
- [[benchmarking-datasets]]
- [[molecular-dynamics-simulation]]
- [[docking-failure-modes]]
- [[protein-ligand-docking]]
- [[benchmarking-datasets]]
- [[virtual-screening]]
- [[autodock-vina]]
- [[molecular-docking-algorithms]]
- [[morehead-2026-evaluating-dl-docking]]
- [[li-2025-decoding-dl-docking-limits]]
- [[van-pinxteren-2025-ml-fep-workflows]]
