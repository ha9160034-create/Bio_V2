# Structure-Based Drug Design

**Summary**: Structure-based drug design (SBDD) uses three-dimensional protein structures to find, optimize, and rationalize small molecules that bind a therapeutic target, with docking, molecular dynamics, and free-energy calculations as its main computational tools.
**Sources**: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Protein–Ligand Docking in the Machine-Learning Era.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulations and drug discovery.md; Basic docking — Autodock Vina 1.2.0 documentation.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Molecular dynamics simulation for all.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md
**Last updated**: 2026-11-01

---

## Definition and place in drug discovery

Structure-based drug design is the use of atomic-resolution structural information about a drug target to guide the discovery and optimization of ligands. Docking is its central fast tool: it predicts how a small molecule binds and provides a score that can be used to rank compounds (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

SBDD is one part of a longer pipeline. A new medicine is estimated to require 5,000–10,000 compounds in chemical and biological studies, about 250 in preclinical testing, and 5 in clinical trials, with 10–15 years and US$1.5–2.0 billion to reach the market (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). Computational structure-based methods aim to reduce the number of compounds that must be made and tested.

## The structural starting point

SBDD requires a three-dimensional receptor structure. Sources include (source: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md):

- X-ray crystallography, historically preferred at resolutions better than 2.5 Å, though modern MD relaxation allows accurate modeling from lower-resolution structures (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md);
- cryo-electron microscopy and NMR;
- homology modelling when no experimental structure exists, with 30%–50% identity models reaching crystal-level accuracy after MD equilibration (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md);
- deep-learning structure prediction, encompassing monomer tools (AlphaFold2, RoseTTAFold) as well as modern all-atom biomolecular co-folding foundation models (AlphaFold 3, Chai-1, Boltz-1) that predict holo complexes directly (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

A co-crystal structure of the target with a ligand is especially valuable because it defines the binding site and provides a validation reference for redocking (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Canonical SBDD workflow

1. **Target selection and validation.** Choose a protein whose modulation should alter disease biology. Most targets are cytosolic kinases, proteases, or nuclear receptors, or membrane-embedded GPCRs and ion channels (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
2. **Structure acquisition and preparation.** Remove unwanted waters, cofactors, and ions; add hydrogens; assign protonation states; and check altlocs and missing residues (source: Basic docking — Autodock Vina 1.2.0 documentation.md; Protein–Ligand Docking in the Machine-Learning Era.md).
3. **Binding-site identification.** Use the co-crystal ligand, a known substrate site, or a pocket predictor; consider allosteric and cryptic sites (source: Protein–Ligand Docking in the Machine-Learning Era.md).
4. **Ligand preparation.** Generate 3D conformers from SMILES or SDF, enumerate protonation and tautomer states, and assign charges appropriate to the docking program (source: Protein–Ligand Docking in the Machine-Learning Era.md).
5. **Docking and scoring.** Perform flexible-ligand/rigid-receptor docking, often across an ensemble of receptor conformations, then rescore or filter (source: Protein–Ligand Docking in the Machine-Learning Era.md).
6. **Experimental testing.** Order or synthesize the top-ranked compounds, typically up to about 100, and assay them (source: Protein–Ligand Docking in the Machine-Learning Era.md).
7. **Lead optimization.** Use structure-guided design, free-energy calculations, and kinetics prediction to improve potency, selectivity, and residence time (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## What docking contributes

Docking supports three SBDD tasks (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md):

- **binding-mode prediction** — where and how a ligand binds;
- **affinity prediction** — ranking analogues during lead optimization;
- **virtual screening** — selecting hits from large compound libraries.

It is also used for structure–activity relationship interpretation, mutagenesis studies, and scaffold hopping (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

The central limitation is accuracy. In the CELPP blind exercise, sampling succeeded 80.4% of the time but scoring only 46.3%, meaning the correct pose was often generated but not ranked first (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). Protein flexibility is identified as the largest remaining hurdle (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## Beyond docking: dynamics and free energy

Static docking treats the receptor as a single structure, but proteins are dynamic ensembles. Molecular dynamics contributes to SBDD by (source: Molecular dynamics simulations and drug discovery.md; Molecular dynamics simulation for all.md):

- identifying **cryptic and allosteric pockets** not visible in crystal structures;
- generating receptor conformations for **ensemble docking** and the **relaxed complex scheme**;
- testing and refining ligand poses;
- predicting **binding free energies** by FEP, TI, MM-PBSA, or MM-GBSA;
- predicting **binding and unbinding kinetics**, including residence time.

The accuracy hierarchy runs from docking scoring functions (fast, approximate) through MM-PBSA/MM-GBSA (intermediate) to alchemical FEP/TI (slow, often ~1 kcal/mol accuracy for congeneric series) (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[binding-free-energy]].

## Success stories

- **Imatinib for chronic myelogenous leukemia.** Imatinib inhibits the Bcr-Abl tyrosine kinase and exemplifies rational structure-based design; it is the standard redocking example in the AutoDock Vina tutorial (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Basic docking — Autodock Vina 1.2.0 documentation.md).
- **Raltegravir for HIV integrase.** MD simulations revealed a cryptic trench on HIV integrase that was absent from crystal structures. Crystallography later confirmed that known inhibitors bind there, and further development at Merck produced raltegravir, the first FDA-approved HIV integrase inhibitor (source: Molecular dynamics simulations and drug discovery.md).
- **FEP-guided HIV reverse transcriptase inhibitors.** A 5 µM non-nucleoside reverse transcriptase inhibitor was optimized to a 55 pM compound using Monte Carlo FEP with OPLS and TIP4P water (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Challenges

- **Scoring accuracy** is target- and test-set-dependent; no universal scoring function exists (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).
- **Protein flexibility** requires ensemble or soft-docking approximations that increase cost and false positives (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Water and entropy** are treated crudely, yet both can dominate binding (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Protonation and tautomer states** can change poses and scores (source: Basic docking — Autodock Vina 1.2.0 documentation.md).
- **Prospective validation is rare.** Most simulation methods are applied retrospectively; standardized benchmarks, blind competitions, and common force fields are still needed (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Polypharmacology.** Drugs often act on several targets, so single-target SBDD can miss efficacy or toxicity mechanisms (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Related pages
- [[protein-ligand-docking]]
- [[virtual-screening]]
- [[binding-site-detection]]
- [[cryptic-and-allosteric-pockets]]
- [[molecular-dynamics-simulation]]
- [[binding-free-energy]]
- [[autodock-vina]]
- [[ma-2024-docking-strategies]]
- [[decherchi-2020-binding-thermodynamics-kinetics]]
