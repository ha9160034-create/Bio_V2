# Ma, Ajibade & Zou (2024) — Docking Strategies

**Summary**: A review of protein–ligand docking that covers the physical basis of molecular recognition, conformational sampling and scoring, docking programs, structure-based and inverse virtual screening, and the field's persistent accuracy problems.
**Sources**: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md
**Last updated**: 2026-11-01

---

## Citation and scope

Zhiwei Ma, Abeeb Ajibade, and Xiaoqin Zou, "Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design," *Communications in Information and Systems* 24(3):199–230 (2024), doi:10.4310/cis.241021221101. Author manuscript PMC11583305, PMID 39584017 (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

The review is written from the perspective of the Zou group, developers of MDock and ITScore. It targets an audience interested in both the physics of recognition and the practical use of docking in CADD.

## Core argument

Docking is framed as solving a three-dimensional jigsaw puzzle: find the best association state between two molecules and estimate how tightly they bind. At a more stringent level, the task is to predict the accurate bound state from the atomic coordinates of the unbound partners (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

The review's central claim is that docking has become indispensable but remains limited: sampling is now much better than scoring, and protein flexibility is the biggest remaining hurdle (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## What the review covers

### 1. Physical basis

- The four main non-covalent interaction classes: hydrogen bonds, ionic interactions, van der Waals forces, and hydrophobic interactions.
- Hydrogen bonds are about 5 kcal/mol, versus about 110 kcal/mol for an O–H covalent bond; van der Waals contacts are about 1 kcal/mol.
- Enthalpy–entropy compensation: strong interactions lower enthalpy but restrict motion, so ΔG changes less than either term.
- Molecular recognition models: lock-and-key, induced fit, conformational selection, and an extended model combining selection with subsequent rearrangement.

(source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md)

### 2. Sampling

A ligand has three translational, three rotational, and one torsional degree of freedom per rotatable bond; proteins add far more. The review describes three ligand-sampling strategies (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md):

- **incremental construction** — build the ligand fragment by fragment in the pocket;
- **multiple conformer generation** — generate low-energy conformers first, then dock each;
- **stochastic search** — Monte Carlo/simulated annealing or genetic algorithms.

Protein flexibility is handled indirectly by ensemble docking, soft docking, or flexible side chains.

### 3. Scoring

The review groups scoring functions into force-field-based, empirical, knowledge-based, and machine-learning categories, and discusses target-specific and consensus scoring. It emphasizes that no universal scoring function exists and that most are sensitive to their training set (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

### 4. Docking programs

More than 60 docking programs have been published. Most treat the ligand as flexible and the protein as rigid. Examples discussed include SLIDE, X-CSCORE, ConsDock, GalaxyDock (flexible side chains), MDock (ensemble docking), ADAM (enlarged pocket), SCARE (alanine scanning to enlarge the pocket), ITScore (retrained on 2,897 complexes), and Vina-carb (carbohydrate terms; 74% versus 55% docking accuracy) (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

### 5. Virtual screening and inverse virtual screening

Conventional VS searches a compound library for binders to one target. Inverse VS (IVS) docks one ligand against many protein targets to find off-targets, alternative indications, or repurposing opportunities. The MDock IVS web server contains 3,268 protein–ligand complexes (3,349 binding sites, 537 targets); a typical run takes 1–4 hours, and a progesterone-receptor example took about one hour on 24 cores (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

The review notes that structure-based and ligand-based VS are complementary, not competitive; a ligand-based pre-screen can reduce the library before docking (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

### 6. Challenges

The most quoted result is from CELPP: **sampling succeeded 80.4% of the time, scoring only 46.3%**. The authors attribute scoring errors to solvent effects, entropic effects, protonation and tautomer states, and other neglected physics. Protein flexibility is the biggest hurdle because sampling all relevant conformations is computationally immense (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## Key numbers

| Quantity | Value |
|---|---|
| Hydrogen-bond strength | ~5 kcal/mol |
| O–H covalent bond | ~110 kcal/mol |
| van der Waals contact | ~1 kcal/mol |
| CELPP sampling success | 80.4% |
| CELPP scoring success | 46.3% |
| ITScore retraining set | 2,897 complexes |
| Vina-carb accuracy | 74% vs 55% for Vina |
| MDock IVS server | 3,268 complexes, 3,349 sites, 537 targets |
| Typical IVS run time | 1–4 hours |

(source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md)

## Strengths

- Connects physics to practice: the interaction and thermodynamics sections explain *why* scoring is hard.
- Covers inverse docking and drug repurposing, which many docking reviews neglect.
- Gives concrete program examples and performance numbers.
- Explicitly states the sampling-versus-scoring gap rather than presenting only success stories.

## Limitations

- The review is perspective-like, not a systematic benchmark; numbers come from different studies and protocols.
- It focuses on the authors' own tools (MDock, ITScore) for several examples.
- Machine learning is covered relatively briefly compared with dedicated reviews.
- No new experimental data.

## Why it matters for this wiki

This is the main source for [[protein-ligand-docking]], [[molecular-docking-algorithms]], [[scoring-functions]], [[protein-ligand-interactions]], [[molecular-recognition-models]], [[virtual-screening]], and [[structure-based-drug-design]]. The CELPP 80.4%/46.3% result is the wiki's central evidence that scoring, not sampling, is the bottleneck.

## Related pages
- [[protein-ligand-docking]]
- [[molecular-docking-algorithms]]
- [[scoring-functions]]
- [[virtual-screening]]
- [[protein-ligand-interactions]]
- [[molecular-recognition-models]]
- [[structure-based-drug-design]]
