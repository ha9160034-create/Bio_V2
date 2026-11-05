# Virtual Screening

**Summary**: Virtual screening ranks large compound libraries computationally to nominate a small set of candidates for experimental testing; structure-based virtual screening does this by docking each molecule into a target binding site and scoring the resulting complex.
**Sources**: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Basic docking — Autodock Vina 1.2.0 documentation.md; Molecular dynamics simulations and drug discovery.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; The Light and Dark Sides of Virtual Screening.pdf; Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf; Application of MM-PBSA Methods in Virtual Screening.pdf; Decoding the limits of deep learning in molecular docking for drug discovery.pdf; Integrating Machine Learning into Free Energy Perturbation Workflows.pdf
**Last updated**: 2026-11-01

---

## Definition and scope

Virtual screening (VS) is the computational search of large chemical libraries for structures predicted to have a desired property — in drug discovery, the ability to bind a target protein (source: Protein–Ligand Docking in the Machine-Learning Era.md).

There are two complementary branches (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md):

- **Structure-based virtual screening (SBVS)** docks each library member into a three-dimensional receptor structure and ranks the complexes with a scoring function.
- **Ligand-based virtual screening (LBVS)** uses the structural and chemical information of known actives — similarity search, pharmacophores, quantitative structure–activity relationships — without docking.

The two are not competitors. A ligand-based pre-screen is often used to reduce the library before structure-based docking, and parallel selection can outperform simple rank fusion (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

The reverse direction is **inverse virtual screening (IVS)**: one ligand is docked against a database of thousands of protein targets to identify plausible targets, off-target effects, or new uses for approved drugs (drug repositioning). The MDock IVS web server contains 3,268 protein–ligand complexes covering 3,349 binding sites and 537 targets; a typical run takes 1–4 hours, and a progesterone-receptor example took about one hour on 24 Intel Xeon cores (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## Why SBVS became practical

Four developments converged to make SBVS at ultra-large scale feasible (source: Protein–Ligand Docking in the Machine-Learning Era.md):

1. Better sampling and scoring methods improved docking, scoring, and screening performance.
2. GPUs and cloud computing dramatically increased throughput.
3. Structural biology (X-ray, NMR, cryo-EM) and deep-learning structure prediction (AlphaFold2, RoseTTAFold, and all-atom co-folding models like AlphaFold 3 and Chai-1) expanded the set of available receptor structures.
4. Make-on-demand chemical space exploded. As of March 2021, the WuXi GalaXi and Enamine REAL Space collections contained 2.1 billion and 17 billion compounds; by June 2022 these had grown to 4.4 billion and 22.7 billion (source: Protein–Ligand Docking in the Machine-Learning Era.md).

With adequate infrastructure it is practical to screen an ultra-large library and select usually up to about 100 compounds for experimental assay (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## The five-step SBVS workflow

### 1. Obtain the receptor structure and the compound library

Receptor structures come from the Protein Data Bank — more than 120,000 unique protein structures as of the review — or from homology modelling, ab initio modelling, AlphaFold2/RoseTTAFold, or biomolecular co-folding models when no experimental structure exists (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Library choice is equally important. Common sources are approved drugs, natural products, synthesized or purchasable compounds and fragments, DrugBank, and the Human Metabolome Database for repurposing. ZINC contains more than 750 million purchasable compounds, including over 230 million in ready-to-dock 3D formats (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Choosing which structure to dock to is not trivial: results depend on conformation, apo versus holo state, and structure quality. The **screening performance index** uses five calculated terms describing how well a set of structures performs on known actives; co-crystal structures with large ligands tend to score well and are good prospective choices. The method requires labelled data, which may not exist for a novel target (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### 2. Detect or choose the binding site

Usually the site is known — the ATP site of a kinase, a protease active site, or a protein–protein interaction interface defined by a co-crystal ligand. When it is unknown, two options exist (source: Protein–Ligand Docking in the Machine-Learning Era.md):

- **blind docking** searches the entire protein surface, at high computational cost;
- **pocket prediction** uses tools such as AlphaSpace, FTMap, MDpocket, Fpocket, and SiteMap to find concave, druggable cavities.

Pocket prediction is also the route to allosteric and cryptic sites, which are important when resistance mutations arise at the orthosteric site (source: Protein–Ligand Docking in the Machine-Learning Era.md). See [[binding-site-detection]] and [[cryptic-and-allosteric-pockets]].

### 3. Prepare the inputs

Protein preparation begins with protonation-state assignment — PROPKA, H++, or SPORES — followed by addition of hydrogens and partial charges, commonly with PDB2PQR. Water molecules and metal ions can be decisive: explicit water-mediated contacts should be analysed, and metal coordination must be modelled for metalloproteins (source: Protein–Ligand Docking in the Machine-Learning Era.md). See [[receptor-preparation]] and [[water-and-solvation]].

Ligands are usually stored as SMILES. RDKit, Open Babel, Omega, or ConfGen generate 3D coordinates; protonation states must be enumerated because they change net charge and atomic charges. AutoDock uses Gasteiger–Marsili charges, whereas AutoDock Vina does not require atomic charges because its scoring terms are charge-independent (source: Protein–Ligand Docking in the Machine-Learning Era.md). See [[ligand-preparation-and-protonation]].

### 4. Dock

The workhorse protocol is **flexible-ligand/rigid-receptor docking**, repeated over an ensemble of receptor conformations (source: Protein–Ligand Docking in the Machine-Learning Era.md). In standard virtual screening engines such as AutoDock Vina, high throughput is achieved by parallelizing independent docking tasks across CPU cores (`--cpu`), optimizing ligand translational, rotational, and active torsional degrees of freedom within a defined search space (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

Strategies to handle receptor flexibility include:

- **ensemble docking** — dock into multiple crystal, NMR, or MD-derived conformations and keep the best-scoring result. Co-crystal structures with larger ligands tend to give the best ensembles, but too many conformers increase false positives and cost linearly (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- **soft docking** — allow partial steric clashes, for example AutoDockFR, which samples explicitly specified side chains with the AutoDock4 force field plus a genetic algorithm and Solis–Wets local search. It improves pose prediction but increases false positives in VS (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- **consensus docking** — combine poses from several programs; a RMSD-based selection of representative poses from 3–11 docking programs outperformed any single program (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### 5. Rescore and filter

Docking scores are built from simplified terms and are unreliable for affinity ranking (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Protein–Ligand Docking in the Machine-Learning Era.md). As demonstrated across CASF and DUD-E evaluations, scoring functions that predict poses well frequently fail to discriminate true active ligands from property-matched inactives in virtual screening (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md). Post-docking filters improve hit rates at low additional cost (source: Protein–Ligand Docking in the Machine-Learning Era.md):

- **Interaction fingerprints (IFP)** convert the 3D interaction pattern into a bit string; Tanimoto similarity to the co-crystal pose discriminates actives from inactives better than classical scoring.
- **GRIM** (Rescoring by Interaction Graph-Matching) compares docked and reference interaction patterns and performs similarly well.
- **Consensus scoring** intersects the top-ranked molecules from two or three scoring functions; this sharply reduced false positives in p38, IMPDH, and HIV protease case studies. The combined functions should not be highly correlated, or their bias is amplified.
- **Machine-learning rescoring** — AtomNet, vScreenML, ΔVina RF20, ΔVina XGB, SIEVE-Score, RF-Score-VS — outperforms classical scoring on benchmark sets, but there is no guarantee of transfer to novel targets outside the training distribution.
- **Free-energy rescoring** — MM-PBSA/MM-GBSA or alchemical methods on selected poses; more accurate but far more expensive (source: Protein–Ligand Docking in the Machine-Learning Era.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Evaluating screening performance

The goal is to rank true actives above inactives in a library. Standard metrics are (source: Protein–Ligand Docking in the Machine-Learning Era.md):

- **Enrichment factor (EF)** — the accumulated fraction of true binders found in the top α percentile of the ranked database (α = 1%, 5%, 10%) divided by the fraction expected at random (higher EF means better early enrichment) (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **ROC curve and AUC** — true-positive rate versus false-positive rate; AUC = 0.5 is random, 1.0 is perfect, most appropriate when actives and inactives are comparable in number (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Spearman rank correlation** — for affinity ranking within a target series (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Case study: ΔLin_F9 XGB on LIT-PCBA

LIT-PCBA contains 15 targets with experimentally curated actives and inactives from PubChem dose–response assays (source: Protein–Ligand Docking in the Machine-Learning Era.md). The protocol:

1. take PDB structures and define the docking box from the co-crystal ligand (source: Protein–Ligand Docking in the Machine-Learning Era.md);
2. process library SMILES with RDKit to generate and protonate low-energy 3D conformers (source: Protein–Ligand Docking in the Machine-Learning Era.md);
3. dock with Smina and the Lin_F9 scoring function, flexible ligand against rigid receptor (source: Protein–Ligand Docking in the Machine-Learning Era.md);
4. rescore the top five poses with ΔLin_F9 XGB and keep the best (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Results across the 15 targets:

- IFP filtering achieved the best average early enrichment, **EF1% = 7.46** (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- ΔLin_F9 XGB achieved the largest number of targets with **EF1% > 2 (13 of 15)** and an average **EF1% = 5.55**, the best among methods that predict binding affinity (source: Protein–Ligand Docking in the Machine-Learning Era.md);
- template-based methods (FINDSITEcomb2.0, FragSite) and the GNINA CNN model achieved lower enrichment (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Comparisons are complicated by differing docking protocols: Sunseri et al. and Yang et al. obtained different average EF1% with Vina because of differences in conformer generation, box definition, and the number of PDB templates used (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Systematic screening failures and chemical artifacts

A large-scale benchmark of DOCK 3.7 and AutoDock Vina on the full DUD-E dataset revealed key empirical boundaries (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf):
- **Throughput differences**: DOCK 3.7 is on average **23 times faster** than AutoDock Vina, outperforming Vina in early enrichment (EF1%) due to physics-based precomputed grids (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Empirical molecular weight bias**: Vina's empirical scoring function correlates strongly with ligand molecular weight, systematically penalizing low-MW fragments and rewarding large compounds (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Pocket rigidity and loop clashes**: In targets with dynamic loops (e.g., MAPK2 80s loop), rigid docking sterically occludes the pocket; in metalloenzymes (COMT $Mg^{2+}$), omitting explicit partial metal charges results in total active enrichment failure (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf). See [[docking-failure-modes]] and [[benchmark-drug-targets]].
- **Deep learning screening collapse**: When applied to virtual screening, end-to-end deep learning docking models perform poorly in separating true active drugs from property-matched decoys; in extensive benchmarks by Li et al., DL scoring on decoy libraries performed no better than random guessing ($\text{EF1\%} \approx 1.0$), failing to match conventional physics-based scoring (source: Decoding the limits of deep learning in molecular docking for drug discovery.pdf).

### The dark side: false-positive chemical artifacts
Top-scoring computational hits frequently exhibit non-specific activity in experimental validation assays (source: The Light and Dark Sides of Virtual Screening.pdf):
- **PAINS**: Sub-structures that trigger false positives via covalent reaction, redox cycling, or fluorescence interference; filtered prior to screening (source: The Light and Dark Sides of Virtual Screening.pdf).
- **Colloidal aggregators**: Hydrophobic small molecules that form 100–1,000 nm colloidal aggregates, non-specifically sequestering and inhibiting target enzymes. They are identified via Aggregator Advisor and verified experimentally by detergent sensitivity (0.01%–0.1% Triton X-100) and validated through [[experimental-binding-assays]] (ITC/SPR) (source: The Light and Dark Sides of Virtual Screening.pdf).

## Multi-tier rescoring: MM-PBSA and active learning FEP

To eliminate false-positive docking hits without exorbitant computational expenditure, tiered screening cascades are deployed (source: Application of MM-PBSA Methods in Virtual Screening.pdf; Integrating Machine Learning into Free Energy Perturbation Workflows.pdf):
- **Endpoint MM-PBSA/GBSA rescoring**: Rescoring the top 500–2,000 docking hits with 1-trajectory explicit-solvent MD simulations relaxes the complex, filters out unstable poses, and boosts EF1% from ~1.2 to >5.0 (source: Application of MM-PBSA Methods in Virtual Screening.pdf). See [[mm-pbsa-and-mm-gbsa]].
- **Active learning FEP**: Coupling Bayesian surrogate models with alchemical FEP screens multi-million-compound libraries up to 20 times faster than brute-force FEP, recovering potent leads with high accuracy (source: Integrating Machine Learning into Free Energy Perturbation Workflows.pdf). See [[alchemical-free-energy-methods]].

## Why hits are still missed

- **Scoring error is larger than sampling error.** In CELPP, sampling succeeded 80.4% of the time but scoring only 46.3% (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- **Receptor flexibility** is the biggest hurdle; a single rigid structure can discard genuine ligands that bind another conformation (source: Molecular dynamics simulations and drug discovery.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- **Decoy bias** inflates retrospective enrichment. DUD-E decoys are presumed inactive without experimental verification and differ from actives in 2D topology; MUV was designed to avoid analogue bias but is better suited to ligand-based VS (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Benchmark leakage** can make ML scoring functions look better than they are on novel targets (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Prospective validation is rare**, which motivated the CACHE initiative launched in 2022 (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Related pages
- [[protein-ligand-docking]]
- [[structure-based-drug-design]]
- [[docking-failure-modes]]
- [[mm-pbsa-and-mm-gbsa]]
- [[alchemical-free-energy-methods]]
- [[deep-learning-docking]]
- [[molecular-docking-algorithms]]
- [[scoring-functions]]
- [[machine-learning-scoring-functions]]
- [[binding-site-detection]]
- [[cryptic-and-allosteric-pockets]]
- [[ligand-preparation-and-protonation]]
- [[receptor-preparation]]
- [[benchmarking-datasets]]
- [[autodock-vina]]
- [[experimental-binding-assays]]
- [[benchmark-drug-targets]]
- [[gimeno-2019-virtual-screening-pitfalls]]
- [[poli-2020-mm-pbsa-virtual-screening]]
- [[xu-2022-docking-failures-analysis]]
- [[chemical-screening-artifacts]]
