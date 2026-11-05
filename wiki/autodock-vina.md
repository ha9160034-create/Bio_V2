# AutoDock Vina

**Summary**: AutoDock Vina is a fast, open-source docking engine that combines a simple empirical scoring function with a stochastic global search and gradient-based local optimization; version 1.2.0 adds the AutoDock4 and Vinardo force fields, macrocycle and hydrated docking, multiple-ligand docking, batch mode, and Python bindings.
**Sources**: autodock-vina-readthedocs-io-en-latest.pdf; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Protein–Ligand Docking in the Machine-Learning Era.md; Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf
**Last updated**: 2026-11-01

---

## What Vina is

AutoDock Vina is one of the docking engines of the AutoDock Suite, together with AutoDock4 (AD4), AutoDock-GPU, AutoDockFR, and AutoDock-CrankPep. It is described as one of the most widely used and successful docking engines, chiefly because it is easy to use, open source, and up to 100× faster than AD4 (source: autodock-vina-readthedocs-io-en-latest.pdf). However, in ultra-large virtual screens, Vina's on-the-fly local optimization is intermediate in throughput compared to rigid precomputed grid engines (such as DOCK 3.7, which is 23× faster on DUD-E benchmarks) (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).

Vina was designed and implemented by Oleg Trott in the Molecular Graphics Lab and is maintained by the Forli Lab at Scripps Research. Its design philosophy is to avoid requiring the user to understand search internals, cluster results, or manipulate quaternions: the user supplies the structures and the search box, and grid maps and atom charges are not required when using the Vina or Vinardo force fields (source: autodock-vina-readthedocs-io-en-latest.pdf).

The original method publication is Trott and Olson (2010); version 1.2.0 is described by Eberhardt, Santos-Martins, Tillack, and Forli (2021) (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Scoring function

Vina uses a united-atom empirical scoring function: only heavy atoms contribute to the score. The positions of hydrogens in the output are therefore arbitrary, but the hydrogens in the input still determine which atoms are hydrogen-bond donors or acceptors, so correct protonation remains essential (source: autodock-vina-readthedocs-io-en-latest.pdf).

Vina ignores user-supplied partial charges. Electrostatic interactions are handled implicitly through the hydrophobic and hydrogen-bonding terms rather than through an explicit Coulomb term (source: autodock-vina-readthedocs-io-en-latest.pdf). This is why no charge assignment step is needed, in contrast to AutoDock4, which uses Gasteiger–Marsili charges (source: Protein–Ligand Docking in the Machine-Learning Era.md).

The scoring-function weights can be changed in the configuration file. For example, `--weight_hydrogen -1.2` doubles the strength of all hydrogen bonds (source: autodock-vina-readthedocs-io-en-latest.pdf). In the classification used by Ma et al., Vina is an empirical scoring function whose score sums inter- and intramolecular energetic contributions (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Version 1.2.0 offers three force fields selected with `--scoring`:

- `vina` — the default Vina force field, which computes its grids internally;
- `ad4` — the AutoDock4 force field, which requires precomputed AutoGrid4 affinity maps;
- `vinardo` — a re-parameterized Vina-like function with modified interaction terms and atomic radii (source: autodock-vina-readthedocs-io-en-latest.pdf; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

> **Important caveat**: scores from the AutoDock and Vina force fields are **not comparable with each other**. The same complex can score −14.72 kcal/mol with AD4 and −13.23 kcal/mol with Vina (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Search algorithm

A Vina docking calculation consists of a number of independent runs, each starting from a random conformation. Each run performs sequential steps of random perturbation, local optimization using the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm, and Metropolis-style acceptance. Local optimization evaluates the scoring function and its derivatives in position–orientation–torsion coordinates (source: autodock-vina-readthedocs-io-en-latest.pdf).

The number of steps per run is chosen heuristically from ligand size and flexibility and from the number of flexible side chains. The number of independent runs is set by **`exhaustiveness`**; because runs execute in parallel, exhaustiveness also limits parallelism. Each run can remember several promising intermediate results, which are merged, refined, clustered, and sorted into the final output (source: autodock-vina-readthedocs-io-en-latest.pdf).

Consequences of this stochastic design:

- docking is **non-deterministic**; identical inputs can give different poses unless the same random seed is supplied;
- the number of output modes is limited both by `--num_modes` and by `energy_range`;
- the search may fail to find the global minimum even when the scoring function is correct (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Search box

Box dimensions are specified in **ångströms**, not in 0.375 Å grid points as in AutoDock4. Confusing the two is a common user error; Vina warns when the box exceeds 27,000 Å³. The recommended practice is to make the box as small as possible without excluding any part of the binding site, and to avoid boxes larger than 30 × 30 × 30 Å unless exhaustiveness is increased (source: autodock-vina-readthedocs-io-en-latest.pdf).

## The basic docking protocol (imatinib–c-Abl)

The reference tutorial docks imatinib (Gleevec, PDB ligand STI) back into the kinase domain of mouse c-Abl (PDB 1IEP), a target for chronic myelogenous leukemia (source: autodock-vina-readthedocs-io-en-latest.pdf).

1. **Receptor preparation.** `mk_prepare_receptor.py` (Meeko) creates a PDBQT with polar hydrogens and partial charges. The tutorial uses a 20 × 20 × 20 Å box centered at (15.190, 53.903, 16.917). Waters, ligands, cofactors, and ions not needed for docking should be removed; `--delete_residues` automates this (source: autodock-vina-readthedocs-io-en-latest.pdf).
2. **Ligand preparation.** `mk_prepare_ligand.py` converts an SDF to PDBQT. PDB format is strongly discouraged for small molecules because it lacks bond connectivity. Protonation and tautomer states must be checked; Molscrub's `scrub.py` can add hydrogens, generate 3D conformers, and enumerate protonation/tautomer states from a SMILES string (source: autodock-vina-readthedocs-io-en-latest.pdf).
3. **Optional AD4 maps.** For `--scoring ad4`, generate a GPF with `mk_prepare_receptor.py -g` and run `autogrid4`, producing `.maps.fld`, per-atom-type maps, a desolvation map `d.map`, and an electrostatic map `e.map` (source: autodock-vina-readthedocs-io-en-latest.pdf).
4. **Docking.** Imatinib is a challenging ligand; the default exhaustiveness is **8**, and the tutorial recommends **32** for consistent results:

```
vina --receptor 1iep_receptor.pdbqt --ligand 1iep_ligand.pdbqt \
     --config 1iep_receptor.box.txt --exhaustiveness=32 \
     --out 1iep_ligand_vina_out.pdbqt
```

5. **Export.** `mk_export.py` converts PDBQT to SDF with correct bond orders and formal charges, using the SMILES stored in the PDBQT header (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Expected results

With exhaustiveness 32, the top AD4 pose is −14.72 kcal/mol and the top Vina pose is −13.23 kcal/mol; both reproduce the crystallographic pose. At the default exhaustiveness, several poses flipped end to end with less favorable energies may be reported (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Flexible receptor docking

Vina can treat selected receptor side chains as flexible. The receptor is split into a rigid PDBQT and a flexible-side-chain PDBQT. Meeko's `-f A:315` flag selects residue 315 of chain A, and `-a` ignores partially resolved residues (source: autodock-vina-readthedocs-io-en-latest.pdf).

The tutorial cross-docks imatinib into c-Abl from PDB 1FPU with Thr315 flexible. AD4 gives a top score of −14.2 kcal/mol; Vina gives −11.63 kcal/mol (source: autodock-vina-readthedocs-io-en-latest.pdf).

Lack of receptor flexibility is described in the Vina documentation itself as "arguably the greatest limitation" of this class of docking methods (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Other protocols in version 1.2.0

### Multiple-ligand docking

Vina can dock several ligands simultaneously, which is relevant to fragment-based design. The tutorial uses PDE in complex with two stereoisomeric inhibitors (PDB 5X72): only the R-isomer binds one region, while both R and S bind a second site. AD4 gives a top score near −18 kcal/mol and needs the top two pose sets for good overlap; Vina gives −19.04 kcal/mol with excellent overlap in the top pose (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Batch mode

`--batch` can be repeated or used with wildcards to dock many ligands sequentially, and `--dir` sets the output directory. This is distinct from simultaneous multiple-ligand docking. There is no native option for iterating over multiple receptors; a shell or PowerShell loop is required (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Zinc metalloproteins

AutoDock4Zn extends the AD4 force field with a potential describing zinc-coordinating ligands, calibrated on 292 crystal complexes. Tetrahedral zinc pseudo-atoms (type TZ) are placed around the zinc, and a special GPF is generated with `prepare_gpf4zn.py` and `AD4Zn.dat`. The tutorial target is PDB 1S63, with a top AD4 score of −13.5 kcal/mol (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Hydrated docking

Explicit bridging waters can be modeled by decorating the ligand with dummy water atoms, adding a water affinity map (`W.map`), docking with AD4, and post-processing with `dry.py` to classify waters as STRONG, WEAK, or DISPLACED. The method improves fragment docking but is **not suitable for virtual screening** because energies must be normalized across ligands (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Macrocycle docking

Rings are normally rigid in Vina. For macrocycles, one ring bond is broken, dummy atoms are added at the broken ends, and a linear attractive potential restores the bond during docking, sampling ring conformations on the fly (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Python API

The `vina` Python package exposes receptor and ligand loading, map computation, scoring, local minimization, docking, and pose writing. `compute_vina_maps()` computes maps for all 22 Vina atom types when no ligand has been loaded, which is convenient for batch screening (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Failure modes and diagnostics

The documentation lists common reasons for a wrong pose (source: autodock-vina-readthedocs-io-en-latest.pdf):

- box specified in grid points instead of ångströms;
- incorrect ligand or receptor protonation;
- stochastic bad luck (try a different seed);
- a correct scoring minimum that the search cannot find (increase exhaustiveness or shrink the box);
- a scoring minimum that simply is not the crystal pose (repeating will not help);
- poor experimental receptor structure, large induced-fit effects, wrong ring conformations, or a 2D ligand input.

In a systematic large-scale benchmark across 102 targets of the DUD-E dataset, Xu et al. (2022) identified four critical algorithmic vulnerabilities in Vina (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf):
1. **Molecular weight (MW) scoring bias**: Vina's score correlates strongly with MW because additive atom-based contact terms systematically assign lower (more favorable) scores to large molecules, severely penalizing low-MW fragments in virtual screening (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
2. **Ligand torsional strain**: Lacking an explicit internal strain penalty, Vina produced 311 types of irregular, distorted torsions across DUD-E poses (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
3. **Metalloenzyme collapse**: Standard Vina preparation (ADT) removes formal metal charges, causing complete failure on metalloenzymes like COMT unless specialized zinc/metal coordination parameters (e.g., AutoDock4Zn) are used (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf).
4. **Throughput**: DOCK 3.7 was on average **23 times faster** than Vina on DUD-E (source: Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf). See [[docking-failure-modes]].

## Related pages
- [[protein-ligand-docking]]
- [[molecular-docking-algorithms]]
- [[scoring-functions]]
- [[docking-failure-modes]]
- [[deep-learning-docking]]
- [[ligand-preparation-and-protonation]]
- [[receptor-preparation]]
- [[virtual-screening]]
- [[water-and-solvation]]
- [[benchmark-drug-targets]]
- [[vina-docs-basic-docking]]
- [[xu-2022-docking-failures-analysis]]
