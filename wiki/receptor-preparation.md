# Receptor Preparation

**Summary**: Receptor preparation selects a biologically relevant protein structure, repairs missing atoms and alternate conformations, assigns protonation states and charges, decides which waters, ions, and cofactors to keep, and converts the result into the format required by the docking program.
**Sources**: Basic docking — Autodock Vina 1.2.0 documentation.md; Protein–Ligand Docking in the Machine-Learning Era.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulation for all.md; autodock-vina-readthedocs-io-en-latest.pdf; Molecular dynamics simulations and drug discovery.md; Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf; gromacs_protein_ligand_complex_tutorial_lemkul.md; Best Practices for Foundations in Molecular Simulations [Article.pdf
**Last updated**: 2026-11-02

---

## Step 1: Choose the structure

The quality of the starting structure limits everything downstream. Selection criteria (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Protein–Ligand Docking in the Machine-Learning Era.md):

- **Resolution** — preferably better than 2.5 Å for quantitative work (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Holo versus apo** — a co-crystal structure defines the binding site and provides a validation ligand, with larger-ligand structures performing well in ensemble selection (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Completeness** — missing loops, side chains, or terminal residues reduce accuracy, especially near the site (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Biological assembly** — the functional oligomer may differ from the asymmetric unit (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Conformation** — the active, inactive, or intermediate state must match the biology being studied (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Species and mutations** — sequence differences from the human target can change pocket shape (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Benchmark Quality Filters** — in standardized affinity benchmarks (PLBenchmarks), Hahn et al. require high-resolution X-ray structures (≤ 2.0–2.5 Å) with R_free < 0.25 and verified Real-Space Correlation Coefficients (RSCC); electron density maps (2Fo - Fc) must be inspected to ensure binding-site sidechains and ligands are unambiguously defined, and crystal packing contacts must be audited to rule out unphysiological lattice constraints (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). See [[simulation-and-benchmark-best-practices]] and [[hahn-2022-binding-affinity-benchmarks]].

When no experimental structure exists, homology modelling, threading, ab initio prediction, or deep learning models (AlphaFold2, RoseTTAFold, or all-atom co-folding models such as AlphaFold 3, Chai-1, and Boltz-1) can supply a model, though side-chain packing and loop placement carry additional uncertainty (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Step 2: Clean the structure

Remove or decide about:

- **Waters** — remove bulk and non-conserved crystallographic solvent for standard rigid docking; preserve explicitly validated structural bridging waters (e.g., Wat301 in HIV-1 protease) that form key bridging hydrogen-bond networks, or employ hydrated docking protocols with explicit water affinity maps (source: Basic docking — Autodock Vina 1.2.0 documentation.md; autodock-vina-readthedocs-io-en-latest.pdf). See [[water-and-solvation]] and [[benchmark-drug-targets]].
- **Ligands and cofactors** — remove the crystallographic ligand unless it is part of the intended site; keep essential cofactors such as NAD, FAD, heme, or metal ions.
- **Ions** — keep structural and catalytic metals; remove buffer ions and cryoprotectants.
- **Altlocs** — choose the highest-occupancy conformation or the one consistent with the binding site.
- **Hydrogens** — usually absent in X-ray structures and must be added.
- **Missing atoms** — rebuild side chains and loops with a modelling tool; flag regions that cannot be rebuilt.
- **N- and C-termini** — assign charged or neutral states consistent with the construct.

The Vina tutorial uses `mk_prepare_receptor.py --delete_residues` to remove unwanted residues automatically (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

## Step 3: Assign protonation states

Protonation changes hydrogen bonding, charge, and often the preferred ligand pose. Tools include PROPKA, H++, and SPORES (source: Protein–Ligand Docking in the Machine-Learning Era.md).

Key decisions:

- **Histidine** — neutral (HID/HIE), doubly protonated (HIP), or doubly deprotonated. Check hydrogen-bond partners and local environment.
- **Cysteine** — thiol versus thiolate; disulfides are usually oxidized.
- **Aspartate and glutamate** — usually deprotonated at pH 7, but buried carboxylates can be protonated.
- **Lysine** — usually protonated, but buried lysines can be neutral.
- **Active-site residues** — catalytic acids and bases often have shifted pKa values; check the mechanism.
- **Termini** — charged unless modified.

For a large virtual screen, a single protonation state per receptor is standard. For mechanistic work, test alternative states and compare.

## Step 4: Add hydrogens and charges

PDB2PQR is a popular tool for adding hydrogens and assigning partial charges (source: Protein–Ligand Docking in the Machine-Learning Era.md). The charge model must match the docking program:

- **AutoDock4** — Gasteiger–Marsili charges, PDBQT atom types.
- **AutoDock Vina** — no charges required; Vina ignores user-supplied partial charges (source: Basic docking — Autodock Vina 1.2.0 documentation.md).
- **Force-field-based docking** — AMBER, CHARMM, or OPLS charges consistent with the force field.

## Step 5: Convert to the docking format

### AutoDock Vina

`mk_prepare_receptor.py` (Meeko) produces a PDBQT with polar hydrogens and charges:

```
mk_prepare_receptor.py -i receptorH.pdb -o receptor -p
```

For flexible side chains, the receptor is split into a rigid PDBQT and a flexible-side-chain PDBQT:

```
mk_prepare_receptor.py -i 1fpu_receptorH.pdb -o 1fpu_receptor -p -v -g \
  --box_size 20 20 20 --box_center 15.190 53.903 16.917 \
  -f A:315 -a
```

`-f A:315` selects residue 315 of chain A as flexible; `-a` ignores partially resolved residues (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

### AutoDock4 maps

For `--scoring ad4`, generate a GPF with `-g` and run `autogrid4`, producing `.maps.fld`, per-atom-type maps, a desolvation map (`d.map`), and an electrostatic map (`e.map`) (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

### Other programs and MD simulation preparation

- **Glide** uses the protein preparation wizard (Maestro).
- **GOLD** accepts MOL2 or PDB with explicit hydrogens.
- **Molecular Dynamics (MD)** requires building a fully parameterized topology, adding missing loops, assigning force-field atom types, adding neutralizing counter-ions (e.g., $\text{Na}^+$, $\text{Cl}^-$), and solvating in explicit water boxes using setup tools like `tleap`, `pdb2gmx`, or CHARMM-GUI (source: Molecular dynamics simulation for all.md). MD relaxation is also frequently employed to relieve crystal-packing strains that distort binding pockets in raw crystallographic structures before docking (source: Molecular dynamics simulation for all.md).

## Step 6: Define the search box

Box dimensions are in **ångströms** for Vina, not in 0.375 Å grid points as in AutoDock4. Confusing the two is a common error. Vina warns when the box exceeds 27,000 Å³. The recommended practice is the smallest box that contains the entire binding site without cutting off any part of the ligand; avoid boxes larger than 30 × 30 × 30 Å unless exhaustiveness is increased (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Special cases

### Membrane proteins

GPCRs, ion channels, and transporters require a lipid bilayer, appropriate membrane orientation, and often a bound agonist or antagonist to stabilize the active state. Docking without the membrane can distort the site.

### Metalloproteins

Zinc, magnesium, calcium, and iron sites need explicit metal ions and, for zinc, specialized atom types. AutoDock4Zn uses tetrahedral zinc pseudo-atoms (type TZ) and a potential calibrated on 292 crystal complexes; the tutorial target is PDB 1S63 with a top AD4 score of −13.5 kcal/mol (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Flexible receptor docking

Selected side chains can be treated as flexible, but every additional flexible residue increases the search space. The c-Abl tutorial makes Thr315 flexible to model the gatekeeper residue; AD4 gives a top score of −14.2 kcal/mol and Vina −11.63 kcal/mol for the cross-docking case (source: autodock-vina-readthedocs-io-en-latest.pdf).

Lack of receptor flexibility is described in the Vina documentation itself as "arguably the greatest limitation" of this class of docking methods (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Ensemble preparation

For ensemble docking, each receptor conformation must be prepared consistently: same protonation rules, same cofactors, same box. Inconsistent preparation creates artificial score differences. The relaxed complex scheme uses multiple MD-derived conformations and ranks ligands by a statistic of the score spectrum (source: Molecular dynamics simulations and drug discovery.md).

## Validation

Before docking a library (source: autodock-vina-readthedocs-io-en-latest.pdf):

1. **Redock the co-crystal ligand** and check RMSD ≤ 2 Å (source: autodock-vina-readthedocs-io-en-latest.pdf).
2. **Check the box** contains the whole site and the ligand (source: autodock-vina-readthedocs-io-en-latest.pdf).
3. **Visually inspect** the site for clashes, missing side chains, and misplaced waters (source: autodock-vina-readthedocs-io-en-latest.pdf).
4. **Test protonation variants** of key residues if the pose is wrong (source: autodock-vina-readthedocs-io-en-latest.pdf).
5. **Run a small decoy set** to check enrichment before committing to a large screen (source: autodock-vina-readthedocs-io-en-latest.pdf).

## Common failure modes

- Box specified in grid points instead of ångströms.
- Incorrect protonation of ligand or receptor.
- Missing cofactor or metal.
- Water molecule blocking the site.
- Alternate conformation inconsistent with the ligand.
- Crystal packing artifact distorting the pocket.
- Poor experimental structure or large induced-fit effects.
- 2D ligand input with no 3D geometry.

(source: autodock-vina-readthedocs-io-en-latest.pdf)

## Related pages
- [[ligand-preparation-and-protonation]]
- [[autodock-vina]]
- [[virtual-screening]]
- [[binding-site-detection]]
- [[water-and-solvation]]
- [[molecular-dynamics-simulation]]
- [[vina-docs-basic-docking]]
