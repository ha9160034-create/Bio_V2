# Ligand Preparation and Protonation

**Summary**: Ligand preparation converts a stored chemical representation — usually SMILES or SDF — into a three-dimensional, protonated, charge-assigned, docking-ready file with correct stereochemistry, tautomers, and bond orders; errors here silently corrupt every downstream pose and score.
**Sources**: Basic docking — Autodock Vina 1.2.0 documentation.md; Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; autodock-vina-readthedocs-io-en-latest.pdf; gromacs_protein_ligand_complex_tutorial_lemkul.md; The Open Force Field Initiative Open Software and Open Science.pdf; Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf
**Last updated**: 2026-11-02

---

## Why preparation matters

Docking programs need explicit 3D coordinates, hydrogens, and connectivity. A ligand file with wrong protonation, missing hydrogens, or a flat 2D conformation will produce a wrong pose no matter how good the search algorithm is. The Vina documentation lists incorrect ligand protonation among the most common reasons for a wrong pose (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

Because Vina uses a united-atom score in which only heavy atoms contribute, the positions of hydrogens in the output are arbitrary — but the hydrogens in the **input** still determine which atoms are donors and acceptors. Correct protonation therefore remains essential (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

## Input formats

Most virtual libraries store compounds as **SMILES** — a line notation encoding atoms and bonds. SMILES has no 3D coordinates and may omit stereochemistry unless specified. **SDF** stores atoms, bonds, and coordinates, and is the preferred exchange format for 3D structures (source: Protein–Ligand Docking in the Machine-Learning Era.md).

PDB format is strongly discouraged for small molecules because it lacks bond connectivity; PDBQT is the AutoDock format that adds atom types and partial charges (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

## Step 1: Parse and standardize

- Remove salts, counterions, and solvent fragments.
- Choose the relevant tautomer and protonation state.
- Assign stereochemistry explicitly; undefined stereocenters generate mixtures that may dock differently.
- Normalize functional groups and aromaticity.
- Check for duplicates and reactive or aggregating compounds.

## Step 2: Generate 3D conformers

Tools include RDKit, Open Babel, Omega, and ConfGen (source: Protein–Ligand Docking in the Machine-Learning Era.md). A typical protocol:

1. Add hydrogens at the chosen pH.
2. Generate multiple low-energy conformers with force-field minimization (MMFF94 or similar).
3. Retain a representative set or the lowest-energy conformer, depending on the docking program.
4. Check that the conformer has sensible bond lengths, angles, and ring geometry.

For the LIT-PCBA case study, library SMILES were processed with RDKit to generate and protonate low-energy 3D conformers before docking (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Step 3: Protonation and tautomers

Ionization changes net charge, hydrogen-bond donor/acceptor patterns, and molecular shape. A carboxylate and a carboxylic acid, or a neutral and a protonated amine, can bind differently and score differently.

- **Rule-based assignment** at a target pH (often 7.4) is standard for large libraries.
- **pKa prediction** with tools such as PROPKA (for proteins) or dedicated ligand pKa predictors improves accuracy for unusual scaffolds.
- **Tautomer enumeration** is important for heterocycles, amides, and keto–enol systems. The Vina documentation recommends Molscrub's `scrub.py`, which can add hydrogens, generate 3D conformers, and enumerate protonation and tautomer states from a SMILES string (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

For a small series, it is often worth docking the two or three most plausible states and comparing scores. For a million-compound library, this is impractical; a single dominant state per compound is the usual compromise.

## Step 4: Charges

Different programs require different charge models (source: Protein–Ligand Docking in the Machine-Learning Era.md):

- **AutoDock4** uses Gasteiger–Marsili atomic charges.
- **AutoDock Vina** does not require atomic charges; its scoring terms are charge-independent and user-supplied partial charges are ignored.
- **AMBER/CHARMM/OPLS** MD force fields require consistent charges — AM1-BCC or RESP — and ligand parameters from GAFF, CGenFF, or OPLS.

Mixing charge models between ligand and receptor is a common source of error in force-field-based scoring.

## Step 5: Convert to the docking format

### AutoDock Vina / AutoDock4

`mk_prepare_ligand.py` (Meeko) converts an SDF to PDBQT:

```
mk_prepare_ligand.py -i ligand.sdf -o ligand.pdbqt
```

The PDBQT stores atom types (including AutoDock-specific types such as A, C, N, OA, NA, HD, and Zn), partial charges, and the SMILES string in the header. The SMILES is used later by `mk_export.py` to restore bond orders and formal charges when converting back to SDF (source: Basic docking — Autodock Vina 1.2.0 documentation.md).

### Other programs

- **Glide** uses its own ligand preparation and accepts SDF/Maestro files.
- **GOLD** uses MOL2 or SDF with explicit hydrogens.
- **Open Babel** can convert between most formats but may not preserve protonation or bond orders reliably.

## Step 6: Validate

Before docking a library, check a sample:

1. Molecular weight and formula are as expected (source: Basic docking — Autodock Vina 1.2.0 documentation.md).
2. Formal charge and protonation are plausible at the assay pH (source: Basic docking — Autodock Vina 1.2.0 documentation.md).
3. Stereocenters are defined (source: Basic docking — Autodock Vina 1.2.0 documentation.md).
4. Bond lengths and angles are reasonable (source: autodock-vina-readthedocs-io-en-latest.pdf).
5. Aromatic rings are planar and kekulized correctly (source: Basic docking — Autodock Vina 1.2.0 documentation.md).
6. The ligand fits the intended binding site without severe strain (source: Basic docking — Autodock Vina 1.2.0 documentation.md).
7. Redocking a known co-crystal ligand reproduces the crystal pose (RMSD ≤ 2 Å) (source: Basic docking — Autodock Vina 1.2.0 documentation.md; autodock-vina-readthedocs-io-en-latest.pdf).

## Special cases

### Macrocycles

Rings are normally rigid in Vina, but macrocycles can be docked flexibly with a specialized protocol: one ring bond is broken, dummy atoms are added at the broken ends, and a linear attractive potential restores the bond during docking, sampling ring conformations on the fly. Since Meeko v0.3.0, macrocycles are made flexible by default (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Covalent inhibitors

Standard docking does not model covalent bond formation. Either dock the non-covalent recognition element and add the covalent link manually, or use a program with explicit covalent docking. The warhead's reactivity and the target's nucleophile must be modelled explicitly.

### Metal-coordinating ligands

Zinc-binding groups need special atom types and maps. AutoDock4Zn places tetrahedral zinc pseudo-atoms (type TZ) and uses a potential calibrated on 292 crystal complexes (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Fragments

Fragments are small and weakly binding; scoring functions calibrated on drug-like ligands may rank them poorly. Hydrated docking can improve fragment results but is not suitable for screening (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Peptides and macrocyclic peptides

These have many rotatable bonds; standard docking may not sample them adequately. Specialized protocols and programs are usually required.

## Small-molecule parameterization for molecular dynamics

When moving from rigid docking to molecular dynamics simulations or alchemical binding free energy calculations, ligand preparation requires generating an explicit all-atom topology with partial charges and force-field parameters (source: gromacs_protein_ligand_complex_tutorial_lemkul.md):

- **CHARMM General Force Field (CGenFF)**: Explicit hydrogens are added, Sybyl .mol2 files must feature unified residue indices (1 JZ4) and sorted bonds via sort_mol2_bonds.pl, and parameter quality is evaluated via penalty scores: penalties $< 10$ denote reliable parameters, 	ext{--}50$ warrant basic QM validation, and $> 50$ mandate manual QM torsional scans (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). See [[gromacs-protein-ligand-protocol]].
- **AMBER / GAFF / GAFF2**: Ligands are typed via Antechamber with AM1-BCC partial charges and converted to GROMACS .itp format using cpype (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Direct Chemical Perception (OpenFF Sage 2.2.0)**: Bypasses discrete atom typing entirely using SMIRNOFF SMARTS/SMIRKS graph queries, eliminating atom-type explosion and automating topology generation for GROMACS, AMBER, and OpenMM (source: The Open Force Field Initiative Open Software and Open Science.pdf). See [[open-force-field-and-ligand-parameterization]] and [[wang-2024-open-force-field-initiative]].
- **Continuous GNN Parameterization (Garnet)**: Uses message-passing neural networks to predict continuous parameters and MBIS partial charges with double exponential nonbonded potentials, eliminating legacy parameter borrowing (source: Training a force field for proteins and small molecules from.pdf). See [[blanco-gonzalez-2024-garnet-force-field]].
- **Affinity Benchmark Curation**: In free energy benchmarks (PLBenchmarks), Hahn et al. require confirming ligand tautomeric and protonation states at the exact experimental assay pH and designing perturbation graphs around the maximum common substructure (MCS) to minimize alchemical phase-space sampling errors (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). See [[simulation-and-benchmark-best-practices]] and [[hahn-2022-binding-affinity-benchmarks]].

## Library filtering

Before docking, it is common to filter by (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Protein–Ligand Docking in the Machine-Learning Era.md):

- molecular weight, logP, hydrogen-bond donors/acceptors, and rotatable bonds (Lipinski-type rules);
- reactive or unstable groups;
- pan-assay interference compounds (PAINS);
- aggregation-prone or highly lipophilic compounds;
- duplicates and already-known actives.

Filtering reduces false positives and computational cost, but it can also remove genuine hits; the filters should match the target and assay.

## Related pages
- [[receptor-preparation]]
- [[autodock-vina]]
- [[virtual-screening]]
- [[water-and-solvation]]
- [[force-fields]]
- [[open-force-field-and-ligand-parameterization]]
- [[wang-2024-open-force-field-initiative]]
- [[blanco-gonzalez-2024-garnet-force-field]]
- [[gromacs-protein-ligand-protocol]]
- [[simulation-and-benchmark-best-practices]]
- [[hahn-2022-binding-affinity-benchmarks]]
- [[lemkul-gromacs-protein-ligand-tutorial]]
- [[vina-docs-basic-docking]]
