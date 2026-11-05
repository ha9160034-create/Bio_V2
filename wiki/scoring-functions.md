# Scoring Functions

**Summary**: A scoring function estimates the binding tightness of a protein–ligand complex; it is the component that ranks docking poses and is the main source of error in docking-based virtual screening.
**Sources**: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md
**Last updated**: 2026-11-01

---

## What a scoring function must do

Scoring functions serve three related but distinct applications (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md):

1. **Binding-mode identification** — rank the experimentally observed pose above hundreds of thousands of decoy poses.
2. **Binding-affinity prediction** — estimate how tightly a ligand binds, chiefly for lead optimization.
3. **Virtual database screening** — rank known binders above inactive compounds in a large library.

An ideal scoring function would perform well on all three, but most existing functions perform well on only one or two (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

## The four families

| Family | Basis | Examples |
|---|---|---|
| Force-field-based | Physical atomic interactions: Lennard-Jones van der Waals + Coulomb electrostatics, plus bonded terms | DOCK, AutoDock, GOLD, SYBYL/D-Score, SYBYL/G-Score |
| Empirical | Weighted sum of interaction terms fitted to experimental binding affinities | FlexX, Glide, ICM, LUDI, PLP, ChemScore, SCORE, X-Score, Surflex, LigScore |
| Knowledge-based | Statistical pair potentials derived from structural databases via the inverse Boltzmann relation | ITScore, PMF, DrugScore, DFIRE, SMoG, BLEEP, MScore, KScore |
| Machine-learning | Nonlinear models trained on structures and affinities | RF-Score, ΔVina RF20, ΔVina XGB, ET-score, ECIF-GBT, AtomNet, P2Rank-style CNNs, GNNs |

Table compiled from (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## Force-field-based scoring

Force-field scores compute Coulombic electrostatics and Lennard-Jones van der Waals terms for non-covalent interactions; more elaborate forms add hydrogen bonds, solvation, and entropic terms (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

The DOCK scoring function uses Amber-derived parameters and a **distance-dependent dielectric constant** to mimic solvent screening. This is efficient, but it cannot describe **desolvation**: charged groups prefer water while non-polar groups prefer non-aqueous environments. Ignoring desolvation biases the score toward highly charged ligands (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

More rigorous solvent treatments include explicit-water free-energy perturbation (FEP) and thermodynamic integration (TI), which are too expensive for database screening. Continuum alternatives are the **Poisson–Boltzmann/surface area (PB/SA)** and **generalized-Born/surface area (GB/SA)** models, often used for post-scoring. Zou and colleagues used GB/SA with solvent-screened electrostatics and a hydrophobic SASA term, optimizing the weights against experimental affinities; the pairwise GB implementation in SDOCK took about 0.5 s per orientation with minimization on a Silicon Graphics Octane R12000 workstation (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

A fundamental difficulty is that the individual energy terms have different scales and are computed by unrelated methods, so empirical weighting coefficients are required. No universal set of weights exists, and individual free-energy terms may not be additive (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

## Empirical scoring

Empirical functions estimate binding affinity as a weighted sum of terms such as van der Waals energy, electrostatics, hydrogen bonds, desolvation, entropy, hydrophobicity, and rotatable bonds. The coefficients are fitted to protein–ligand complexes with known structures and affinities (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

- **SCORE1** (Böhm) used 45 complexes and four terms: hydrogen bonds, ionic interactions, lipophilic contact surface, and number of rotatable bonds. It was later expanded to 82 complexes and additional terms including binding-site electrostatics, water molecules, cavities, and aromatic interactions.
- **ChemScore** (Eldridge et al.) considered hydrogen bonds, metal atoms, lipophilic effects, and rotatable bonds, calibrated on 82 complexes and tested on sets of 20 and 10.
- **X-Score** (Wang et al.) used 200 complexes and four terms: van der Waals, hydrogen bonds, hydrophobic effects, and effective rotatable bonds.

The main weakness is dependence on the training set; double counting of interactions is a further problem (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

## Knowledge-based scoring

Knowledge-based functions derive pair potentials from the observed frequency of atom pairs in experimental structures using the inverse Boltzmann relation:

ΔG ∝ −k_B T ln [ρ(r)/ρ*(r)]

where ρ(r) is the pair density at distance r and ρ*(r) is the pair density in a reference state with zero interactions (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

These potentials are mean-force potentials rather than true interaction potentials, and the **reference state is the central unsolved problem**. Most implementations use an atom-randomized reference state that ignores excluded volume and connectivity (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

Representative functions:

- **DrugScore** used 17 atom types and 1,376 complexes; a CSD-based version, DrugScore^CSD, followed.
- **BLEEP** used 40 atom types; it reached R = 0.71 on nine serine protease–inhibitor complexes and R = 0.74 on 90 complexes, but only R² = 0.28 on Muegge and Martin's set.
- **PMF** used 34 ligand atom types and 16 protein atom types from 697 structures, with a ligand-volume correction; it reached R² = 0.61 on 77 complexes. PMF04 used 7,152 complexes.
- **DFIRE** reached R = 0.63 on 100 protein–ligand complexes.
- **ITScore** avoids the reference-state problem by iteratively adjusting pair potentials until they reproduce the observed pair distribution functions, including both native structures and decoys. ITScore/SE adds solvation and configurational entropy (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

Knowledge-based functions are fast and relatively insensitive to the training set, but their pairwise potentials are derived from well-bound structures and may score wrong poses well (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

## Consensus scoring

Consensus scoring combines multiple scoring functions to reduce the bias of any single one. Strategies include vote-by-number, number-by-number, rank-by-number, average rank, and linear combination. Examples are MultiScore, X-CScore, GFscore, SCS, and SeleX-CS (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

Consensus works best when the component functions are **not highly correlated**; combining correlated functions amplifies rather than cancels bias (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## How scoring functions are evaluated

### Binding-mode prediction: self-docking versus prospective blind scoring (CELPP)

The standard metric for pose prediction is root-mean-square deviation (RMSD) between the predicted and experimental ligand coordinates; $\text{RMSD} \le 2.0\text{ \AA}$ denotes success (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md). However, reported accuracy depends critically on whether the benchmark is retrospective self-docking or blind prospective scoring (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Protein–Ligand Docking in the Machine-Learning Era.md):
- **Self-docking / Redocking**: In retrospective redocking into native holo co-crystal structures with pre-formed cavities, scoring functions achieve high nominal success rates (typically 70% to 85% of complexes with $\text{RMSD} \le 2.0\text{ \AA}$) (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- **Prospective blind challenges (CELPP)**: In the Continuous Evaluation of Ligand Pose Prediction (CELPP) challenge, scoring functions were evaluated prospectively on newly released PDB complexes (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). CELPP demonstrated that while conformational search algorithms generated near-native poses ($\text{RMSD} \le 2.0\text{ \AA}$) in **80.4%** of cases, scoring functions correctly ranked the near-native pose first in only **46.3%** of targets (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

RMSD can also be confounded by small or symmetric ligands (which can obtain low RMSD by chance) and by large flexible ligands (where solvent-exposed groups inflate RMSD despite a correct binding core) (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). Alternatives include interaction-based classification (IBAC) and real-space R-factor (RSR) (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

### Affinity prediction

The Pearson correlation between predicted scores and experimental binding energies is standard; because the relationship need not be linear, the Spearman rank correlation is also used. As detailed in biophysical studies, experimental uncertainty across assays, buffer ionic strengths, and pH conditions sets an upper limit on achievable correlation (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Virtual screening

Enrichment measures the fraction of known actives recovered above a percentile of the ranked database. ROC/AUC is preferred when the number of inactives is comparable to the number of actives (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md). See [[benchmarking-datasets]].

## Benchmark evidence: no universal winner

On Wang et al.'s 100-complex set, the success rates (RMSD ≤ 2.0 Å) and affinity correlations were (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md):

| Scoring function | Type | Pose success (%) | Affinity R |
|---|---|---|---|
| ITScore/SE | K | 91 | 0.65 |
| DrugScore^CSD | K | 87 | 0.62 |
| ITScore | K | 82 | 0.65 |
| Cerius2/PLP | E | 76 | 0.56 |
| SYBYL/F-Score | E | 74 | 0.30 |
| Cerius2/LigScore | E | 74 | 0.35 |
| DrugScore^PDB | K | 72 | 0.60 |
| Cerius2/LUDI | E | 67 | 0.36 |
| X-Score | E | 66 | 0.64 |
| AutoDock | F | 62 | 0.05 |
| DFIRE | K | 58 | 0.63 |
| DOCK/FF | F | 58 | 0.40 |
| Cerius2/PMF | K | 52 | 0.40 |
| SYBYL/G-Score | F | 42 | 0.56 |
| SYBYL/ChemScore | E | 35 | 0.47 |
| SYBYL/D-Score | F | 26 | 0.48 |

The table shows that pose success and affinity correlation are not the same thing. SYBYL/F-Score achieves 74% pose success but only R = 0.30 for affinity; AutoDock achieves 62% pose success but R = 0.05 on this set (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

Enrichment is similarly inconsistent. At the top 5% of four screened databases, ITScore recovered 19.2% (ERα), 68.3% (MMP3), 34.9% (fXa), and 37.0% (AChE) of actives, against theoretical maxima of 39.2%, 88.3%, 43.7%, and 97.5%. PMF-Score, despite R² = 0.61 in affinity prediction, recovered 0% on ERα and AChE and only 5.0% on MMP3 (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

> **Contradiction to note**: A 2010 review reported that the machine-learning function RF-Score reached R = 0.953 on 1,105 complexes (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md). Later work found that such machine-learning functions can reproduce strong scoring power but may perform poorly on docking and screening (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). The two statements are not necessarily incompatible — they measure different tasks — but they show why a single high correlation cannot establish scoring-function quality.

## Recent directions

- **Target-specific scoring functions**: training on a particular protein family can improve pose and affinity prediction, but such functions do not generalize to novel targets (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- **Added interaction terms**: halogen-bond parameters (AutoDock VinaXB), carbohydrate CH–π terms (Vina-Carb, SLICK), and metal terms.
- **Machine learning**: see [[machine-learning-scoring-functions]].
- **Better benchmarks**: see [[benchmarking-datasets]].

## Persistent challenges

Current scoring functions remain far from universally accurate. Key problems are the treatment of desolvation and entropy, the lack of universal weighting coefficients, transition-metal parameterization, test-set dependence, and the difficulty of comparing functions evaluated on different datasets (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

## Related pages
- [[protein-ligand-docking]]
- [[machine-learning-scoring-functions]]
- [[autodock-vina]]
- [[binding-free-energy]]
- [[water-and-solvation]]
- [[benchmarking-datasets]]
- [[virtual-screening]]
- [[benchmark-drug-targets]]
- [[huang-2010-scoring-functions]]
