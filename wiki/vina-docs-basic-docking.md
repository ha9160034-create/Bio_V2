# AutoDock Vina 1.2.0 Documentation — Basic Docking and Tutorials

**Summary**: The official AutoDock Vina 1.2.0 manual and tutorial set, covering installation, the basic imatinib–c-Abl redocking protocol, flexible receptor docking, multiple-ligand and batch docking, zinc metalloproteins, hydrated docking, macrocycles, and the Python API.
**Sources**: Basic docking — Autodock Vina 1.2.0 documentation.md; autodock-vina-readthedocs-io-en-latest.pdf
**Last updated**: 2026-11-01

---

## What this source is

The AutoDock Vina documentation for release 1.2.0, maintained by the Forli Lab at The Scripps Research Institute. Vina was originally designed and implemented by Oleg Trott in the Molecular Graphics Lab; the code is open source under the Apache License 2.0. The manual describes Vina as a turnkey docking program built on a simple scoring function and rapid gradient-optimization conformational search, requiring only the structures and the search space — no grid maps or atom charges when using the Vina or Vinardo force fields (source: autodock-vina-readthedocs-io-en-latest.pdf).

Vina is described as up to **100× faster** than AutoDock4 and as improving the average accuracy of binding-mode prediction relative to AD4. Derivatives and extensions include QuickVina2 (search), Smina (interface and scoring terms), Vina-Carb (carbohydrates), VinaXB (halogen bonds), and Vinardo (re-parameterized scoring) (source: autodock-vina-readthedocs-io-en-latest.pdf).

Version 1.2.x added functionality previously available only in AutoDock4: macrocycle flexibility, specialized metal-coordination models, explicit-water hydrated docking, simultaneous docking of multiple ligands, and Python bindings (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Structure of the documentation

- Availability and installation (pip, conda, precompiled binaries, source build).
- Participating: FAQ, citations, change log, software requirements.
- Tutorials: basic docking, flexible docking, multiple ligands, batch mode, zinc metalloproteins, hydrated docking, macrocycles, Python scripting.
- Python API reference.

## The basic docking protocol

The tutorial docks imatinib (Gleevec; PDB ligand from entry **1IEP**) into the kinase domain of mouse c-Abl, a target for chronic myelogenous leukemia (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

### 1. Receptor preparation

`mk_prepare_receptor.py` (Meeko) produces a PDBQT with polar hydrogens and partial charges:

```
mk_prepare_receptor.py -i 1iep_receptorH.pdb -o 1iep_receptor -p -v \
  --box_size 20 20 20 --box_center 15.190 53.903 16.917
```

The `-p` flag writes the receptor PDBQT; `-v` with the box arguments writes a config TXT and a box PDB for visualization. The `--delete_residues` option removes unwanted waters, ligands, cofactors, and ions. The ADFR `prepare_receptor` tool is an alternative, but it requires a coordinate file that already contains all hydrogens, which can be added with REDUCE (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

### 2. Ligand preparation

`mk_prepare_ligand.py` converts an SDF to PDBQT. PDB format is strongly discouraged for small molecules because it lacks bond connectivity. If the ligand has no hydrogens, `scrub.py` from Molscrub can protonate it, generate 3D conformers, and enumerate tautomers and protonation states (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

The warning is explicit: always check protonation before docking — "your success can sometimes hang by just an hydrogen atom" (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

### 3. Optional affinity maps for AutoDock4 scoring

To use the AD4 force field, generate a GPF with `-g` and run `autogrid4`. This produces `.maps.fld`, per-atom-type affinity maps (A, C, HD, H, NA, N, OA), a desolvation map `d.map`, and an electrostatic map `e.map`. The Vina and Vinardo force fields do not need these maps (source: Basic docking — Autodock Vina 1.2.0 documentation.md; autodock-vina-readthedocs-io-en-latest.pdf).

### 4. Running the docking

The default exhaustiveness is **8**; the tutorial recommends **32** because imatinib is a challenging ligand and Vina may occasionally miss the correct pose at default settings. AD4 and Vina force fields give different scores, and the documentation warns that energies from different force fields are not comparable (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

Results for the same 1IEP system (source: Basic docking — Autodock Vina 1.2.0 documentation.md):

| Force field | Top score (kcal/mol) |
|---|---|
| AutoDock4 | −14.72 |
| Vina | −13.23 |

The tutorial notes that with exhaustiveness 32 Vina usually returns a single pose at this energy, whereas at the default exhaustiveness several end-to-end flipped poses with less favorable energy may appear (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

### 5. Export

`mk_export.py` converts the PDBQT output back to SDF, restoring bond orders and formal charges from the SMILES stored in the PDBQT header (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

## How the search works

The documentation's FAQ gives an unusually clear description (source: autodock-vina-readthedocs-io-en-latest.pdf):

- Each docking calculation consists of a number of **independent runs starting from random conformations**.
- Each run performs sequential steps: a random perturbation, a local optimization with the **Broyden–Fletcher–Goldfarb–Shanno (BFGS)** algorithm, and an accept/reject decision.
- The number of steps is set heuristically from ligand size and flexibility.
- The number of runs is set by **exhaustiveness**, which therefore also limits parallelism.
- Unlike AutoDock4, each Vina run can retain several promising intermediate results; these are merged, refined, clustered, and sorted automatically.

## Search box

- Box dimensions are in **ångströms**, not in 0.375 Å grid points as in AutoDock4. Confusing the two is a common error (source: autodock-vina-readthedocs-io-en-latest.pdf).
- Use the **smallest box that still contains the whole site and ligand**. Vina warns above **27,000 Å³**, roughly a 30 × 30 × 30 Å box, and larger boxes require higher exhaustiveness (source: autodock-vina-readthedocs-io-en-latest.pdf).
- Ligand and flexible side-chain atoms outside the box are not sampled.

## Common causes of a wrong pose

The FAQ lists the usual suspects (source: autodock-vina-readthedocs-io-en-latest.pdf):

1. Box specified in grid points instead of ångströms.
2. Incorrect ligand or receptor protonation.
3. Bad luck with the random seed — retry with a different seed.
4. The scoring minimum is correct but the search cannot find it — raise exhaustiveness or shrink the box.
5. The scoring minimum is simply wrong — repeated runs may occasionally give the right answer by accident.
6. Poor-quality X-ray or NMR receptor structure.
7. Large induced-fit effects not captured by the rigid receptor.
8. Rigid rings in the wrong conformation.
9. A flat 2D ligand input.
10. The experimental bound conformation differs from the crystal or NMR model.

Vina uses a **united-atom** scoring function: only heavy atoms contribute, so hydrogen positions in the output are arbitrary. Input hydrogens still determine donor/acceptor assignment, so protonation remains essential (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Flexible receptor docking

The tutorial cross-docks imatinib into c-Abl from PDB **1FPU**, treating **Thr315** (the gatekeeper residue) as flexible. The receptor is split into rigid and flexible PDBQT files:

```
mk_prepare_receptor.py -i 1fpu_receptorH.pdb -o 1fpu_receptor -p -v -g \
  --box_size 20 20 20 --box_center 15.190 53.903 16.917 \
  -f A:315 -a
```

Results (source: Basic docking — Autodock Vina 1.2.0 documentation.md):

| Force field | Top score (kcal/mol) |
|---|---|
| AutoDock4 | −14.2 |
| Vina | −11.63 |

The documentation states plainly that **lack of receptor flexibility is arguably the greatest limitation** of this class of docking methods (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Multiple-ligand docking

Vina 1.2.0 can dock several ligands simultaneously, which is useful for fragments, cofactors, or two inhibitors occupying adjacent sites. The tutorial uses PDE in complex with two inhibitors (PDB **5X72**) (source: autodock-vina-readthedocs-io-en-latest.pdf).

```
vina --ligand 5x72_ligand_p59.pdbqt 5x72_ligand_p69.pdbqt \
  --maps 5x72_receptor --scoring ad4 --exhaustiveness 32 \
  --out 5x72_ligand_ad4_out.pdbqt
```

Results: AD4 top score **−17.67 kcal/mol**; Vina top score **−19.04 kcal/mol**, with the best Vina pose overlapping the crystallographic coordinates for one isomer (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Batch mode

For sequential screening, `--batch` accepts one or more ligand files and can be combined with shell wildcards or a loop over receptors and ligands. This is the command-line route to virtual screening; the Python API is the programmatic alternative (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Zinc metalloproteins

The tutorial uses **1S63** and the **AutoDock4Zn** force field. `zinc_pseudo.py` adds tetrahedral zinc pseudo-atoms (type **TZ**) around the zinc ion, representing preferred coordination positions. `prepare_gpf4zn.py` writes a GPF with the `AD4Zn.dat` parameter file, and `autogrid4` produces the maps (source: autodock-vina-readthedocs-io-en-latest.pdf).

The top AD4 score is **−13.5 kcal/mol** and corresponds to the crystallographic pose. AutoDock4Zn's potential was calibrated on 292 crystal complexes (source: autodock-vina-readthedocs-io-en-latest.pdf; Basic docking — Autodock Vina 1.2.0 documentation.md).

## Hydrated docking

Some waters are so strongly bound and conserved that they behave as part of the receptor. The hydrated-docking protocol decorates the ligand with dummy water atoms, adds a water affinity map (`W.map`), docks with AutoDock4, and post-processes with `dry.py`, which classifies waters as **STRONG**, **WEAK**, or **DISPLACED** (source: autodock-vina-readthedocs-io-en-latest.pdf).

The tutorial docks the fragment nicotine into acetylcholine-binding protein (**1UW6**), where removing waters can cause inaccurate scoring or an incorrect pose. The top AD4 score is **−8.261 kcal/mol** (source: autodock-vina-readthedocs-io-en-latest.pdf).

Two caveats are explicit:

- the method was calibrated with the **AutoDock4** force field, and the documentation advises against using it with Vina or Vinardo;
- energies must be normalized across chemically diverse ligands, so the method is **not suitable for virtual screening**, although it is fine for comparing poses of the same ligand.

## Macrocycles

Vina cannot directly treat ring-bond flexibility. The specialized protocol (source: autodock-vina-readthedocs-io-en-latest.pdf):

1. break one ring bond, producing an open form with independent torsions;
2. add a dummy atom to each broken end, preserving the original bond length and 1–3 angle;
3. apply a linear attractive potential between the dummy atoms during docking to restore the closed ring.

This samples macrocycle conformations on the fly, at the cost of extra rotatable bonds. Since **Meeko v0.3.0**, macrocycles are made flexible by default, so the basic docking tutorial suffices; `mk_export.py` restores correct connectivity and removes dummy atoms from the output (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Python API

The Python bindings expose the docking engine directly (source: autodock-vina-readthedocs-io-en-latest.pdf):

```python
from vina import Vina

v = Vina(sf_name='vina')          # or 'ad4', 'vinardo'
v.set_receptor('1iep_receptor.pdbqt')
v.set_ligand_from_file('1iep_ligand.pdbqt')
v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])

energy = v.score()
energy_minimized = v.optimize()
v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)

v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('1iep_ligand_vina_out.pdbqt', n_poses=5, overwrite=True)
```

Notes:

- `cpu=0` uses all available CPUs.
- Flexible side chains are supplied as a second receptor file: `v.set_receptor('rigid.pdbqt', 'flex.pdbqt')`.
- If no ligand is loaded before `compute_vina_maps()`, maps are computed for all **22** Vina atom types, which is convenient for batch screening of unknown libraries.
- `dock(exhaustiveness=32, n_poses=20)` runs 32 Monte Carlo samplings and stores up to 20 poses; output is limited by `energy_range`, default 3 kcal/mol.

## What the documentation emphasizes

- **Ease of use**: few parameters, no charges or maps for the default force field, automatic clustering and ranking.
- **Speed**: multicore execution and up to 100× faster than AD4.
- **Reproducibility**: the random seed is printed; supplying the same seed reproduces a run exactly, while changing the number of CPUs or exhaustiveness can change results.
- **Honesty about failure**: the FAQ lists scoring errors, protonation, box errors, rigid rings, and induced fit as reasons docking fails.
- **Force-field incomparability**: AD4, Vina, and Vinardo scores are not on the same scale and should not be compared directly.

## Strengths

- Complete, command-level recipes that can be reproduced exactly.
- Covers the practical edge cases — zinc, water, macrocycles, multiple ligands, batch mode — that general docking reviews omit.
- The FAQ explains the algorithm and failure modes better than most papers.
- The Python API makes large-scale and custom workflows accessible.

## Limitations

- Tutorial-centric: it teaches how to run Vina, not how to judge whether a docking result is biologically meaningful.
- The default scoring function is not calibrated for absolute affinity; scores are ranking heuristics.
- No guidance on receptor conformational ensembles or cryptic sites beyond flexible side chains.
- Performance numbers and version-specific behavior change between releases.
- The Python bindings are documented as Linux and Mac only in this release, although the command-line executable runs on Windows.
- **Empirical biases identified in modern benchmarks**: In an independent 102-target DUD-E evaluation, Xu et al. (2022) revealed that Vina's empirical scoring function correlates strongly with molecular weight (penalizing low-MW fragments), generates 311 irregular torsion types due to omitting an explicit internal strain penalty, and collapses on catalytic metalloenzymes (such as COMT $\text{Mg}^{2+}$) unless specialized coordination parameters (such as AutoDock4Zn) are provided. See [[docking-failure-modes]] and [[benchmark-drug-targets]].

## Why it matters for this wiki

This is the primary source for [[autodock-vina]], and a supporting source for [[ligand-preparation-and-protonation]], [[receptor-preparation]], [[binding-site-detection]], [[water-and-solvation]], [[protein-ligand-docking]], and [[virtual-screening]]. It supplies the concrete protocol parameters — box center and size, exhaustiveness, force-field choice, PDBQT preparation, and the 1IEP/1FPU/5X72/1S63/1UW6 examples — that make the docking pages reproducible.

## Related pages
- [[autodock-vina]]
- [[ligand-preparation-and-protonation]]
- [[receptor-preparation]]
- [[binding-site-detection]]
- [[water-and-solvation]]
- [[protein-ligand-docking]]
- [[virtual-screening]]
- [[scoring-functions]]
- [[docking-failure-modes]]
- [[benchmark-drug-targets]]
