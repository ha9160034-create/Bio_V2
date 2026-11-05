# Poziemski et al. (2026) — MD Time-Series Descriptors in Affinity Prediction

**Summary**: A large-scale benchmark evaluating machine learning models trained on molecular dynamics time-series descriptors across 2,752 protein–ligand complexes, demonstrating how dynamic features improve binding affinity prediction and generalizability over static crystal structures.
**Sources**: `Assessment of molecular dynamics time series.pdf`
**Last updated**: 2026-11-02

---

## The dynamic dimension in machine learning scoring

Poziemski et al. address a fundamental limitation of conventional structure-based machine learning (ML) and deep learning (DL) scoring functions: the overwhelming majority of models are trained exclusively on static three-dimensional crystallographic coordinates (source: Assessment of molecular dynamics time series.pdf). Static structures ignore thermal fluctuations, induced-fit conformational plasticity, water residence times, and the entropic persistence of non-covalent interactions (source: Assessment of molecular dynamics time series.pdf).

While molecular dynamics (MD) trajectories contain rich information regarding complex stability and fluctuation dynamics, previous attempts to incorporate MD snapshots into ML scoring often degraded predictive power (source: Assessment of molecular dynamics time series.pdf). The authors identify that incorporating raw, uncurated MD frames introduces severe simulation noise and high-dimensional redundancy, which drowns out the physical signal in ML models (source: Assessment of molecular dynamics time series.pdf).

## The large-scale MD dataset and time-series extraction

To systematically evaluate the value of MD dynamics, the authors constructed the largest curated MD dataset for ML affinity prediction to date:
- **Dataset Scale**: **2,752 protein–ligand complexes** curated from the PDBbind database with experimentally determined binding affinities ($\text{p}K_d / \text{p}K_i$) (source: Assessment of molecular dynamics time series.pdf).
- **Simulation Protocol**: Explicit-solvent all-atom MD trajectories generated under standardized, reproducible equilibration and production protocols (source: Assessment of molecular dynamics time series.pdf).
- **Time-Series Feature Extraction**: Trajectory time series were extracted using the `tsfresh` library across 75 distinct biophysical feature types, including interatomic contact distances, hydrogen-bond counts, solvent-accessible surface areas (SASA), hydrophobic burial, and Lennard-Jones/Coulombic energy components (source: Assessment of molecular dynamics time series.pdf).
- **Extracted Dynamic Properties (`ts_descriptors`)**: Features capture statistical, time-domain, and frequency-domain behaviors, such as mean values, variance, standard deviation, autocorrelation lengths, skewness, kurtosis, fast Fourier transform (FFT) spectral coefficients, and contact persistence fractions (source: Assessment of molecular dynamics time series.pdf).

## Five-stage feature curation and trimming pipeline

To eliminate simulation noise and collinearity, the authors implemented a rigorous five-step dimensionality reduction workflow:
1. **Initial Extraction**: Calculating statistical and time-series features across MD trajectories for all 75 biophysical property types (source: Assessment of molecular dynamics time series.pdf).
2. **$p$-Value Trimming**: For parametric and non-parametric time-series descriptors within each feature type, the variant with the lowest $p$-value against experimental affinity is retained (source: Assessment of molecular dynamics time series.pdf).
3. **Correlation Pruning**: Pairwise Pearson correlation coefficients (PCC) are evaluated; for any descriptor pair with $\text{PCC} \ge 0.8$, the feature with lower cardinality is discarded to prevent multi-collinearity (source: Assessment of molecular dynamics time series.pdf).
4. **Mutual Information Maximization**: Descriptors are filtered by mutual information against experimental binding affinity, selecting at most one optimal `ts_descriptor` per underlying physical metric (source: Assessment of molecular dynamics time series.pdf).
5. **Final Representation**: Each complex is characterized by a compact vector consisting of at most **51 MD `ts_descriptors`**, **19 static ligand property descriptors**, **63 static crystallographic descriptors**, and a **1024-bit ECFP4 topological fingerprint** (source: Assessment of molecular dynamics time series.pdf).

## Experimental validation and generalization regimes

Models—including **Random Forest (RF)**, **XGBoost (XGB)**, and **Multilayer Perceptrons (MLP)**—were evaluated across three distinct data-splitting paradigms (source: Assessment of molecular dynamics time series.pdf):
- **Random Split (80:20)**: Baseline cross-validation where ligands and targets are randomly distributed; demonstrates high baseline performance across all architectures (source: Assessment of molecular dynamics time series.pdf).
- **Target Split (UniProt Generalization)**: Complexes are grouped strictly by UniProt accession ID, ensuring the test set contains only completely unseen protein families; tests the true transferability of learned biophysical features (source: Assessment of molecular dynamics time series.pdf).
- **Scaffold Split (Bemis–Murcko Generalization)**: Ligands are partitioned into clusters based on their Bemis–Murcko core scaffolds using DeepChem, placing less abundant and distinct chemical matter into the test set to evaluate chemotype out-of-distribution performance (source: Assessment of molecular dynamics time series.pdf).

## Key findings and implications for CADD

- **Noise Suppression is Essential**: Feeding raw snapshot descriptors directly into ML models degrades scoring accuracy as trajectory length increases, confirming that unpruned thermal fluctuations act as noise (source: Assessment of molecular dynamics time series.pdf).
- **Superior Generalization on Challenging Splits**: When curated via the five-stage trimming pipeline, MD time-series descriptors significantly outperform purely static crystallographic models on both **Target Splits** and **Scaffold Splits** (source: Assessment of molecular dynamics time series.pdf).
- **Physical Insights**: The most informative `ts_descriptors` relate to contact persistence (how long key protein–ligand hydrogen bonds or salt bridges remain intact over time) and the amplitude of pocket solvent fluctuation, parameters that are entirely inaccessible from single static crystal structures (source: Assessment of molecular dynamics time series.pdf).

## Related pages
- [[machine-learning-scoring-functions]]
- [[molecular-dynamics-simulation]]
- [[scoring-functions]]
- [[benchmarking-datasets]]
- [[binding-free-energy]]
- [[virtual-screening]]
- [[simulation-and-benchmark-best-practices]]
- [[gromacs-protein-ligand-protocol]]

