# Xu et al. (2022) — Docking Failures in Large-Scale Virtual Screening

**Summary**: A large-scale benchmark of DOCK 3.7 and AutoDock Vina on the DUD-E dataset systematically dissecting the causes of docking and virtual screening failures.
**Sources**: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Min Xu, Cheng Shen, Jincai Yang, Qing Wang, and Niu Huang, "Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening," *ACS Omega* 7(44):39417–39428 (2022), doi:10.1021/acsomega.2c05826 (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

The authors systematically analyze the success and failure mechanisms of structure-based virtual screening (SBVS) by executing head-to-head calculations between a grid-based physics-driven docking engine (UCSF DOCK 3.7) and a widely used empirical stochastic engine (AutoDock Vina) across the 102 benchmark targets of the Directory of Useful Decoys, Enhanced (DUD-E) (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

## Quantitative performance comparison

The benchmark establishes significant computational and enrichment distinctions between the two programs (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf):
- **Computational efficiency**: On average, DOCK 3.7 was **23 times faster** than AutoDock Vina across the DUD-E dataset using default operating parameters (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Early enrichment (EF1%)**: DOCK 3.7 outperformed Vina in early recognition (top 1% enrichment factor, EF1%) across the majority of targets, which the authors attribute to precomputed physics-based electrostatic and van der Waals potential grids (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Target family dependency**: Performance varied by target architecture: Vina excelled on nuclear receptors (such as the androgen receptor, ANDR), whereas DOCK 3.7 outperformed on metalloenzymes (such as catechol-O-methyltransferase, COMT) and deep hydrophobic cavities (such as the $\beta_2$-adrenergic receptor, ADRB2) (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- **Active chemotype complementarity**: The top-scoring 1% of actives selected by DOCK 3.7 and Vina showed only a minor fraction of overlapping compounds, indicating that their distinct search and scoring mechanisms capture complementary chemical spaces (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

## The five primary docking failure mechanisms

Through detailed residue-level analysis of six representative targets—ADRB2, MAPK2, FKB1A, ANDR, FA10, and COMT—the authors categorized the root causes of docking failures (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf):

### 1. Molecular weight (MW) scoring bias
Vina's empirical scoring function correlates strongly with molecular weight: because its scoring terms accumulate additively with atom count, high-MW molecules receive disproportionately negative (favorable) kcal/mol scores regardless of binding complementarity (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf). Consequently, low-MW and fragment-like actives are systematically ranked poorly in Vina screens. In contrast, DOCK 3.7's grid-based force-field score shows little correlation with ligand MW, avoiding this size penalty (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

### 2. Ligand internal strain energy and distorted torsions
Neither program successfully enforces rigorous conformational strain penalties during docking search (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf):
- Vina produced 311 types of irregular, energetically unfavorable torsions in its top poses, while DOCK 3.7 generated 281 types (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- Vina poses frequently featured distorted ring geometries and strained rotatable bonds because its empirical objective function lacks an explicit torsional strain energy term (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

### 3. Protein pocket rigidity and loop clashes
The rigid receptor approximation was identified as a dominant failure mode when ligand entry or binding requires protein conformational rearrangement (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf):
- In MAPK2, the mobile 80s loop sterically occluded the ATP-binding pocket in the static crystal structure, causing true actives to clash and receive unfavorable scores in both programs (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- In ANDR, the rigid precalculated conformer ensembles of DOCK 3.7 suffered internal clashes with the narrow, tight receptor cavity, whereas Vina's iterated local search successfully optimized flexible torsions within the confining boundary (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

### 4. Metalloenzyme coordination errors
Metalloenzymes represent a major vulnerability for empirical docking (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf):
- In COMT, which requires coordination with a catalytic $Mg^{2+}$ ion, AutoDock Vina completely failed (near-zero enrichment) because standard AutoDockTools (ADT) preparation removes formal charges and ignores explicit metal coordination parameters (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
- DOCK 3.7 successfully enriched COMT actives by redistributing 0.2 electron units to each coordinating active-site atom.
- When Vina was supplemented with AutoDock4Zn-derived parameters for coordinated metal ions, its adjusted logAUC improved dramatically to 18.29, proving that explicit metal charge modeling is required to prevent screening failure (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

### 5. Multi-subpocket scoring discrepancies
In targets featuring multiple extended subpockets (e.g., coagulation factor Xa, FA10), poses often established favorable interactions in one subpocket while clashing or mispairing in adjacent subpockets, exposing scoring function imbalances in penalizing partial desolvation (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

## Related pages
- [[virtual-screening]]
- [[autodock-vina]]
- [[protein-ligand-docking]]
- [[scoring-functions]]
- [[benchmarking-datasets]]
- [[receptor-preparation]]

