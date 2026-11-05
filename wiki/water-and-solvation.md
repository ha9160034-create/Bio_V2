# Water and Solvation

**Summary**: Water is not a passive background: it determines the hydrophobic effect, mediates protein–ligand contacts, creates desolvation barriers, and contributes entropy and enthalpy that can dominate binding; docking and free-energy methods treat it only approximately.
**Sources**: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Basic docking — Autodock Vina 1.2.0 documentation.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulation for all.md; Protein–Ligand Docking in the Machine-Learning Era.md; autodock-vina-readthedocs-io-en-latest.pdf; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md
**Last updated**: 2026-11-01

---

## Why water matters

A protein–ligand–solvent system is a thermodynamic system of solute (protein and ligand) and solvent (water and buffer ions). The driving forces for association are a synthetic result of interactions and heat exchange among protein, ligand, water, and ions (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Before binding, both protein and ligand are hydrated. Forming a complex requires disrupting protein–water and ligand–water interactions, releasing water, and reorganizing the solvent around the new interface. These processes contribute to both ΔH and ΔS (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## The hydrophobic effect

Nonpolar groups tend to exclude water and aggregate. The classical explanation is an entropy gain from releasing ordered water molecules at the solute surface, although the mechanism remains debated; scaled-particle theory treats the effect as dependent on solute size (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

In the binding entropy decomposition, the solvent term ΔS_solv is usually positive and favorable because surface burial releases constrained water (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). The lock-and-key model is entropy-dominated precisely because perfectly matched rigid surfaces displace a large amount of water without a conformational entropy cost (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Water as part of the interaction network

- **Bridging waters** mediate hydrogen bonds between protein and ligand, extending the contact network and sometimes determining specificity (source: autodock-vina-readthedocs-io-en-latest.pdf; Protein–Ligand Docking in the Machine-Learning Era.md).
- **Crystallographic waters** are often reproduced in MD simulations and can be functionally important. In the Src kinase simulations, crystallographic water positions were consistently recovered, and a water shell around PP1 created an entrance kinetic barrier (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Water-mediated allostery** — changes in hydration can propagate structural effects.
- **Proton wires** — chains of water and titratable residues can transfer protons, relevant to enzymes and membrane channels.

## Desolvation barriers

A ligand entering a buried site must shed its hydration shell. In the β-adrenergic receptor simulations, the highest free-energy barrier was not from the extracellular vestibule to the orthosteric site but from bulk solvent to the vestibule. Analysis showed this barrier was determined by near-complete desolvation of the ligand and partial desolvation of the binding pocket (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Desolvation is therefore a major determinant of both k_on and k_off, and a source of enthalpy–entropy compensation: releasing water is entropically favorable but breaking water–solute interactions is enthalpically costly (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Heat capacity and surface burial

The change in heat capacity ΔC_p correlates strongly with buried surface area. Burial of apolar surface gives negative ΔC_p; burial of polar surface gives positive ΔC_p. Because hydration water has different properties from bulk water, ΔC_p links thermodynamics to structure and reports on desolvation (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## How docking handles water

Most docking protocols remove water entirely, which is a major approximation (source: Basic docking — Autodock Vina 1.2.0 documentation.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

Options in practice:

### Remove all waters

Standard for large-scale virtual screening: fastest, but loses bridging waters and water-mediated contacts.

### Keep selected waters

Retain waters that mediate known interactions or occupy the site. This requires judgment and structural evidence; keeping all waters blocks the site.

### Hydrated docking (AutoDock)

Explicit bridging waters can be modeled by decorating the ligand with dummy water atoms, adding a water affinity map (W.map), docking with AutoDock4, and post-processing with `dry.py`. Waters are classified as STRONG, WEAK, or DISPLACED and scored. The method improves fragment docking but is **not suitable for virtual screening** because energies must be normalized across ligands (source: autodock-vina-readthedocs-io-en-latest.pdf).

### Explicit-water rescoring

MM-PBSA with explicit waters near the ligand can improve correlation with experiment; the improvement depends on the system and the starting pose. In a bromodomain study, including explicit waters changed the Pearson correlation, but not consistently (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Water in machine-learning scoring

ΔVina XGB includes three bridge-water features, and ΔLin_F9 XGB includes bridge-water terms, reflecting the importance of water in affinity prediction (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## How MD handles water

MD normally uses explicit water. The water model — TIP3P, TIP4P, SPC, and variants — affects density, dielectric behavior, diffusion, and the balance of hydrophobic and polar interactions. The water model must be consistent with the force field and ion parameters (source: Molecular dynamics simulation for all.md).

Implicit solvent models (Poisson–Boltzmann, Generalized Born, SASA) are used in MM-PBSA/MM-GBSA to reduce cost. They treat water as a continuum and lose explicit hydrogen bonding and bridging waters. PB is theoretically superior to GB, but studies have found the opposite in some systems, and both depend on radii and charges (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Water entropy is hard to compute

The entropic contribution of water in the binding pocket is often not negligible and is difficult to model. It is a major source of error in MM-PBSA and in free-energy calculations, and it is one reason docking scores correlate poorly with affinity (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

In a host–guest study of β-cyclodextrin, the entropy gain from releasing water absorbed on the cavity was a major driving force for binding, and ΔG estimated from kinetics was more accurate than ΔG from phase-space averages (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Practical checklist

1. Inspect the structure for ordered waters at the site; decide which to keep based on B-factors, hydrogen-bonding role, and conservation.
2. For virtual screening, remove waters but flag sites where bridging waters are known.
3. For pose validation or free-energy work, use explicit solvent and check whether key waters remain stable.
4. When using MM-PBSA/GBSA, test whether including explicit waters changes the ranking.
5. Report the water model and water treatment; results depend on them.

## Related pages
- [[protein-ligand-interactions]]
- [[binding-thermodynamics]]
- [[binding-free-energy]]
- [[molecular-dynamics-simulation]]
- [[force-fields]]
- [[receptor-preparation]]
- [[autodock-vina]]
- [[biochemistry-free-for-all]]
