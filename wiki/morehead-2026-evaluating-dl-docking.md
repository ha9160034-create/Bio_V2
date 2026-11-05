# Morehead et al. (2026) — Assessing Deep Learning for Protein–Ligand Docking

**Summary**: A Nature Machine Intelligence study systematically evaluating deep-learning docking and co-folding models against conventional physics-based docking on physical validity and PoseBusters benchmarks.
**Sources**: Assessing the potential of deep learning for protein–ligand docking.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Alex Morehead, Nabin Giri, Jian Liu, Pawan Neupane, and Jianlin Cheng, "Assessing the potential of deep learning for protein–ligand docking," *Nature Machine Intelligence* 8:32–41 (2026), doi:10.1038/s42256-025-01160-1 (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

The article conducts a comprehensive, multi-benchmark assessment of emerging deep learning (DL) docking models and all-atom biomolecular co-folding systems, measuring their performance against conventional physics-based docking engines on structural accuracy, physical validity, and chemical realism (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

## Methods and benchmark architectures evaluated

The authors benchmarked nine computational pipelines across primary and multiligand docking scenarios (source: Assessing the potential of deep learning for protein–ligand docking.pdf):
- **Conventional physics-based baseline**: P2Rank (pocket detection) paired with AutoDock Vina (v.1.2.5) for molecular docking.
- **Deep learning docking algorithms**: DiffDock-L, DynamicBind, NeuralPLexer, and RoseTTAFold-All-Atom (RFAA).
- **All-atom biomolecular co-folding models**: AlphaFold 3 (AF3), Chai-1, and Boltz-1 (tested with both full Multiple Sequence Alignments [MSAs] and single-sequence inputs).
- **Evaluation datasets**: PoseBusters Benchmark ($n=130/308$ complexes), Astex Diverse Set, DockGen-E (unseen receptor folds), and CASP15 Protein–Ligand Interface (PLI) targets (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

## The PoseBusters validity metric (PB-Valid)

A critical contribution of the study is moving beyond raw root-mean-square deviation (RMSD) to enforce chemical and physical validity using the PoseBusters suite (PB-Valid) (source: Assessing the potential of deep learning for protein–ligand docking.pdf). PB-Valid assesses:
- Steric clashes between protein and ligand atoms;
- Interligand steric clashes in multi-binding pockets;
- Aromatic ring planarity and non-aromatic ring puckering;
- Standard covalent bond length and bond angle deviations;
- Chirality preservation and internal torsional strain (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

The central metric is the **combined success rate**, defined as the percentage of predicted complexes achieving both structural accuracy ($\text{RMSD} \le 2.0\text{ \AA}$) and full chemical validity (PB-Valid) (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

## Key empirical findings

### 1. The "RMSD illusion" and clash rates in generative models
Generative diffusion and regression DL models (such as DiffDock-L and DynamicBind) often report high raw pose accuracy rates based solely on heavy-atom centroid RMSD $\le 2$ Å (source: Assessing the potential of deep learning for protein–ligand docking.pdf). However:
- Their PB-Valid rates are severely depressed because their coordinate generation frequently introduces severe atomic overlaps, unphysical bond lengths, and distorted rings.
- Consequently, their true combined success rate ($\text{RMSD} \le 2.0\text{ \AA}$ and PB-Valid) drops steeply, proving that raw RMSD without geometric validation produces an exaggerated estimate of practical docking efficacy (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

### 2. Co-folding models and the MSA dependency
All-atom co-folding systems (Chai-1, Boltz-1, and AlphaFold 3) demonstrated superior structural performance compared to earlier DL docking tools on primary ligand prediction (source: Assessing the potential of deep learning for protein–ligand docking.pdf):
- However, their predictive power depends heavily on evolutionary sequence information. When evaluated in single-sequence mode (AF3-Single-Seq, Boltz-1-Single-Seq, Chai-1-Single-Seq), performance dropped substantially.
- This demonstrates that co-folding algorithms rely heavily on evolutionary conservation patterns in the MSA to infer binding pocket geometries rather than directly recognizing small-molecule physical forces (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

### 3. Competitiveness of conventional physics-based docking
The conventional baseline—P2Rank paired with AutoDock Vina—consistently exhibited **100% PB-Valid compliance** because it explores rigid bond parameters within empirical force-field energy landscapes (source: Assessing the potential of deep learning for protein–ligand docking.pdf). On the combined success metric ($\text{RMSD} \le 2.0\text{ \AA}$ and PB-Valid), P2Rank-Vina remained highly competitive with or outperformed several specialized DL docking models across the benchmark sets (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

### 4. Multiligand and cofactor vulnerability
In multi-ligand complexes (co-factors, catalytic metal ions, or multiple substrates), generative DL tools exhibited frequent interligand steric clashes, whereas physics-based docking programs model multi-body electrostatic and Lennard-Jones repulsion naturally (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

### 5. Posthoc physical relaxation
Applying posthoc energy minimization to DL-predicted poses resolved atomic clashes but frequently caused the ligand centroid to shift out of the pocket, revealing that the generative model's learned probability density does not align with the physical thermodynamic free energy minimum (source: Assessing the potential of deep learning for protein–ligand docking.pdf).

## Related pages
- [[protein-ligand-docking]]
- [[machine-learning-scoring-functions]]
- [[benchmarking-datasets]]
- [[autodock-vina]]
- [[virtual-screening]]
- [[molecular-docking-algorithms]]

