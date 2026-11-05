# Yang, Chen & Zhang (2022) — Protein–Ligand Docking in the ML Era

**Summary**: A review of machine-learning scoring functions and structure-based virtual screening, covering datasets, evaluation metrics, ML architectures from random forests to graph neural networks, the VS workflow, and the LIT-PCBA benchmark.
**Sources**: Protein–Ligand Docking in the Machine-Learning Era.md
**Last updated**: 2026-11-01

---

## Citation and scope

Chao Yang, Eric Anthony Chen, and Yingkai Zhang, "Protein–Ligand Docking in the Machine-Learning Era," *Molecules* 27(14):4568 (2022), doi:10.3390/molecules27144568. PMC9323102, PMID 35889440 (source: Protein–Ligand Docking in the Machine-Learning Era.md).

The review focuses on **machine-learning scoring functions** and **structure-based virtual screening**. It is the wiki's main source for modern ML scoring, benchmarking datasets, and the LIT-PCBA case study.

## Structure of the review

1. Introduction to docking and scoring.
2. Datasets and evaluation metrics.
3. Machine-learning scoring functions.
4. Structure-based virtual screening, including the docking protocol and workflow.
5. Case study on LIT-PCBA.
6. Challenges.

## Datasets

The review organizes datasets into three groups (source: Protein–Ligand Docking in the Machine-Learning Era.md):

### Structure + affinity

- **PDBbind** — version 2020 general set has **19,443 complexes** with K_d, K_i, or IC50 from **1.2 pM to 10 mM**, updated annually. It also provides a refined set and a core set used by CASF.
- **CSAR** — Community Structure–Activity Resource exercises.
- **D3R Grand Challenge** — blind pose and affinity competitions.

### Active/inactive labels

- **DUD-E** — 102 targets, 22,886 actives, 50 computer-generated decoys per active. Decoys have similar physicochemical properties but dissimilar 2D topology; many are presumed inactive without experimental verification, a major drawback.
- **MUV** — 17 targets, 30 actives and 15,000 experimentally verified inactives per target, designed to avoid analogue bias. Better suited to ligand-based VS because many actives are not validated against the putative target.
- **LIT-PCBA** — 15 targets from PubChem dose–response assays with homogeneous conditions and careful removal of false positives (0.5 < Hill slope < 2.0). Limitation: 8 of 15 targets are cell-based phenotypic assays.

### Bioactivity only

- **BindingDB** — 2,513,948 binding data for 8,839 targets and 1,077,922 small molecules (as of 4 May 2022).
- **ChEMBL** — version 30 has 19,286,751 activities for 14,885 targets and 2,157,379 compounds.

## Evaluation metrics

- **Pose prediction** — RMSD ≤ 2 Å is success; use the minimum symmetry-corrected RMSD.
- **Affinity prediction** — Pearson correlation, or Spearman rank correlation when the relationship is monotonic but not linear.
- **Screening** — enrichment factor (EF) at 1%, 5%, or 10%; ROC curve and AUC. AUC = 0.5 is random, 1.0 is perfect, and AUC is most appropriate when actives and inactives are comparable in number.
- **CASF-2016** — 57 targets × 5 complexes (285 total), evaluated on scoring, ranking, docking, and screening. Its screening assessment is limited by the small database and unverified inactives.

(source: Protein–Ligand Docking in the Machine-Learning Era.md)

## Machine-learning scoring functions

### Early methods

- **RF-Score** (Ballester and Mitchell, 2010) — random forest on protein–ligand atom-type pair counts within a distance cutoff. The first ML scoring function to outperform classical functions on scoring.
- **SFCscore^RF** (2013) — random forest on 63 empirical descriptors (ligand-dependent, interaction, and surface terms).
- Both performed much worse on docking and screening than classical scoring functions.

(source: Protein–Ligand Docking in the Machine-Learning Era.md)

### Delta-machine-learning

Instead of predicting affinity from scratch, Δ-ML learns a correction to a physical score. This preserves the docking power of the base function while improving affinity accuracy (source: Protein–Ligand Docking in the Machine-Learning Era.md).

- **ΔVina RF20** (2017) — 10 pharmacophore SASA features + 10 of Vina's 58 empirical terms; random forest; best on CASF 2007 and 2013 across scoring, ranking, docking, and screening.
- **ΔVina XGB** (2019) — Vina's 58 features + 30 SASA + 3 bridge-water + 2 ligand-stability + 1 metal-count feature; XGBoost; trained on solvated and dry structures plus docking decoys; better on CASF-2016.
- **ΔLin_F9 XGB** (2022) — Gaussian terms characterizing interactions at different distance ranges; enlarged training set with weak binders and docking poses; superior scoring, ranking, and screening on CASF-2016.

### Interaction-feature methods

- **ET-score** (2021) — distance-weighted interatomic contacts between atom-type pairs; Pearson R = 0.827 on CASF-2016.
- **ECIF-GBT** (2021) — extended connectivity interaction features; Pearson R = 0.857 on CASF-2016.
- Both were trained solely on crystal structures, and their docking and screening performance remains an issue.

### Deep learning

- **NNScore 1.0/2.0** (2010/2011) — neural networks; 1.0 classifies active/inactive from 194 features, 2.0 estimates pK_d.
- **AtomNet** (2017) — first CNN scoring function using 3D grids over the interface; four 3D convolutional layers and two fully connected layers; better AUC than Vina on DUD-E.
- **Pafnucy, Kdeep, OnionNet** — CNN variants. Pafnucy and Kdeep use atom-property or atom-type grids; OnionNet uses rotation-free element-pair contacts in shells, addressing the coordinate-frame dependence of earlier CNNs.
- **PotentialNet** (2018) — graph convolutional network with covalent-only, dual noncovalent/covalent, and ligand-gathering stages; Pearson R = 0.822 on CASF-2007.
- **graphDelta, graphBAR, SIGN** — later GNNs that did not consistently outperform traditional ML scoring functions on CASF.

### Target-specific scoring

Models trained on a single target or family can reach state-of-the-art performance on well-studied targets with sufficient data, but cannot be applied to a novel target with little experimental data (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Structure-based virtual screening

### Why it grew

Four developments: better sampling and scoring; GPUs and cloud computing; more structures from X-ray, NMR, cryo-EM, AlphaFold2, and RoseTTAFold; and enormous make-on-demand libraries. WuXi GalaXi grew from 2.1 billion to 4.4 billion compounds and Enamine REAL Space from 17 billion to 22.7 billion between March 2021 and June 2022. ZINC contains more than 750 million purchasable compounds, including over 230 million ready-to-dock 3D structures (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### Five-step workflow

1. Obtain receptor structure and compound library.
2. Detect or choose the binding site — co-crystal, blind docking, or pocket prediction (AlphaSpace, FTMap, MDpocket, Fpocket, SiteMap).
3. Prepare inputs — protonation (PROPKA, H++, SPORES), hydrogens and charges (PDB2PQR), waters and metals, ligand 3D conformers (RDKit, Open Babel, Omega, ConfGen).
4. Dock — flexible ligand/rigid receptor, often over an ensemble; soft docking (AutoDockFR) or multi-program consensus.
5. Rescore and filter — IFP, GRIM, consensus scoring, ML rescoring, free-energy rescoring.

(source: Protein–Ligand Docking in the Machine-Learning Era.md)

### Handling protein flexibility

- **Ensemble docking** — multiple crystal, NMR, or MD conformations; too many conformers increase false positives and cost linearly.
- **Soft docking** — allow partial clashes; AutoDockFR samples explicit side chains with the AutoDock4 force field plus GA and Solis–Wets local search. Better poses but more false positives.
- **MD** — find cryptic pockets, generate ensembles, and validate complexes.

## LIT-PCBA case study

The ΔLin_F9 XGB protocol (source: Protein–Ligand Docking in the Machine-Learning Era.md):

1. Use PDB structures with co-crystal ligands to define the box.
2. Process library SMILES with RDKit into protonated low-energy 3D conformers.
3. Dock with Smina and the Lin_F9 scoring function, flexible ligand against rigid receptor.
4. Rescore the top five poses with ΔLin_F9 XGB; keep the best.

Results:

- **IFP** filtering gave the best average early enrichment, **EF1% = 7.46**.
- **GRIM** performed similarly.
- **ΔLin_F9 XGB** gave **13 of 15 targets with EF1% > 2** and **average EF1% = 5.55**, the best among affinity-predicting methods.
- Template-based methods (FINDSITEcomb2.0, FragSite) and GNINA gave lower enrichment.

The comparison is complicated by different docking protocols — conformer generation, box definition, and number of templates — which changed Vina-based EF1% between studies.

## Challenges identified

- Scoring accuracy remains the bottleneck; no universal function exists.
- Protein flexibility limits both sampling and scoring.
- Decoy bias and false negatives inflate retrospective enrichment.
- ML models may not generalize to novel targets outside the training distribution.
- Most validation is retrospective; CACHE (2022) was launched to provide prospective testing.
- Ultra-large libraries require efficient conformer generation and docking infrastructure.

## Why it matters for this wiki

This is the main source for [[machine-learning-scoring-functions]], [[benchmarking-datasets]], [[virtual-screening]], and much of [[scoring-functions]]. It supplies the wiki's modern dataset numbers and the LIT-PCBA comparison.

## Related pages
- [[machine-learning-scoring-functions]]
- [[benchmarking-datasets]]
- [[virtual-screening]]
- [[scoring-functions]]
- [[protein-ligand-docking]]
- [[ligand-preparation-and-protonation]]
- [[receptor-preparation]]
