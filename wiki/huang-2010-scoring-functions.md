# Huang, Grinter & Zou (2010) — Scoring Functions and Their Evaluation

**Summary**: A perspective reviewing the three classical scoring-function families, consensus scoring, the criteria used to evaluate pose prediction, affinity prediction, and screening, and the databases available for benchmarking.
**Sources**: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Scoring functions and their evaluation methods for protein-ligand docking recent advances and future directions.pdf; Protein–Ligand Docking in the Machine-Learning Era.md
**Last updated**: 2026-11-01

---

## Citation and scope

Sheng-You Huang, Sam Z. Grinter, and Xiaoqin Zou, "Scoring functions and their evaluation methods for protein-ligand docking: recent advances and future directions," *Physical Chemistry Chemical Physics* 12(40):12899–12908 (2010), doi:10.1039/c0cp00151a. Author manuscript PMC11103779, PMID 20730182 (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md; Scoring functions and their evaluation methods for protein-ligand docking recent advances and future directions.pdf).

The paper is a perspective, not a new benchmark. Its value is the systematic framing of what a scoring function must do, how it is built, and how it should be tested.

## The three applications

The review opens by distinguishing three applications that are often conflated (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md):

1. **Binding-mode and site determination** — rank the experimental pose above hundreds of thousands of decoys.
2. **Absolute affinity prediction** — estimate binding tightness, chiefly for lead optimization.
3. **Virtual database screening** — rank known binders above inactives.

An accurate function should perform well on all three, but most perform well on only one or two. This is the paper's organizing insight.

## The three scoring-function families

### Force-field-based

These compute physical atomic interactions: Lennard-Jones van der Waals and Coulomb electrostatics, plus bond stretching, bending, and torsional terms. Parameters come from experimental data and ab initio quantum mechanics (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

The central difficulty is solvent treatment. DOCK uses a distance-dependent dielectric to mimic screening, but this cannot describe desolvation: charged groups prefer water, nonpolar groups prefer non-aqueous environments. Ignoring desolvation biases scores toward highly charged ligands. PB/SA and GB/SA continuum models are more rigorous but expensive; a pairwise GB implementation in SDOCK took about 0.5 s per orientation with minimization on a 1990s SGI Octane R12000 workstation (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

A further problem is that the individual energy terms have different scales and are computed by unrelated methods, so empirical weights are required; no universal weights exist, and the terms may not be additive (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

### Empirical

Empirical functions are weighted sums of terms such as van der Waals, electrostatics, hydrogen bonds, desolvation, entropy, hydrophobicity, and rotatable bonds, fitted to complexes with known affinities (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

Examples:

- **SCORE1** (Böhm) — 45 complexes; hydrogen bonds, ionic interactions, lipophilic contact surface, rotatable bonds. Later expanded to 82 complexes with electrostatics, water, cavities, and aromatic terms.
- **ChemScore** (Eldridge) — hydrogen bonds, metals, lipophilic effects, rotatable bonds; calibrated on 82 complexes.
- **X-Score** (Wang) — 200 complexes; van der Waals, hydrogen bonds, hydrophobic effects, effective rotatable bonds.

Weaknesses: dependence on the training set and double counting of interactions.

### Knowledge-based

Knowledge-based functions derive pair potentials from the observed frequency of atom pairs in experimental structures using the inverse Boltzmann relation:

$$\Delta G \propto -k_B T \ln\frac{\rho(r)}{\rho^*(r)}$$

where ρ(r) is the observed pair density and ρ*(r) is the reference density (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

The potentials are mean-force potentials, not true interaction potentials, and the **reference state is the central unsolved problem**. Most implementations use an atom-randomized reference that ignores excluded volume and connectivity (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

Examples include ITScore, PMF, DrugScore, DFIRE, SMoG, BLEEP, MScore, GOLD/ASP, and KScore.

### Consensus scoring

Consensus scoring combines several functions to reduce the bias of any one. The review cautions that the combined functions should not be highly correlated, or their bias is amplified rather than cancelled (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

## Evaluation criteria

The paper distinguishes (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md):

- **Pose prediction** — RMSD between predicted and experimental pose; success at ≤2 Å, with symmetry correction for equivalent atoms.
- **Affinity prediction** — Pearson correlation with experimental values; Spearman rank correlation when the relationship is monotonic but not linear.
- **Screening** — enrichment factor at a fixed percentage of the ranked database, ROC curves, and AUC.

It also discusses the composition of a good test set: high-quality structures, a wide range of proteins and affinities, drug-like non-covalent ligands, and no training/test overlap.

## Benchmark results discussed

On a set of 100 diverse complexes, performance varied enormously (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md):

- **ITScore/SE** ranked 91% of native poses correctly (RMSD ≤ 2.0 Å).
- **SYBYL/D-Score** ranked only 26%.
- For affinity correlation, **ITScore/SE and ITScore** reached R = 0.65, while **AutoDock** reached R = 0.05 on that set.

These numbers are test-set dependent and should not be read as universal rankings. They illustrate the paper's point that no single function dominates all tasks.

The review also describes an enrichment benchmark with 146 ERα, 60 MMP3, 129 fXa, and 54 AChE inhibitors, each set mixed with 999 random drug-like molecules as inactives (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

## Databases discussed

The paper reviews the protein–ligand databases available at the time for developing and validating scoring functions, including PDBbind and CASF, CSAR, and the DUD decoy sets (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md). See [[benchmarking-datasets]] for the modern versions and their biases.

## Limitations

- Published in 2010; the specific programs and benchmark numbers are historical. Modern ML scoring functions, deep-learning docking, and ultra-large libraries are not covered.
- Retrospective self-docking bias: High pose success rates (e.g. 91% for ITScore/SE) reflect redocking into pre-formed co-crystal cavities. In subsequent prospective blind evaluations organized by the Zou group in the CELPP challenge, prospective scoring succeeded in only **46.3%** of targets (despite conformational sampling succeeding in 80.4%), demonstrating that static retrospective benchmarks significantly overestimate prospective scoring accuracy (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- The benchmark numbers come from the authors' own and others' studies with different protocols, so cross-comparisons are informal.
- The reference-state problem in knowledge-based scoring is stated but not solved.
- The paper predates the DUD-E false-negative critique and the LIT-PCBA benchmark.

## Why it matters for this wiki

This is the foundational source for [[scoring-functions]] and [[benchmarking-datasets]]. Its three-task framing and its warning that a function can be good at one task and bad at another explain why RF-Score can excel at affinity prediction yet fail at docking and screening (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Related pages
- [[scoring-functions]]
- [[machine-learning-scoring-functions]]
- [[benchmarking-datasets]]
- [[protein-ligand-docking]]
- [[yang-2022-ml-docking]]
