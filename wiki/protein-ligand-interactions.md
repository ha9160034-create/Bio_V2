# Protein–Ligand Interactions

**Summary**: Protein–ligand binding is stabilized by non-covalent interactions — hydrogen bonds, ionic pairs, van der Waals contacts, hydrophobic burial, and weaker effects such as CH–π and halogen bonds — whose net enthalpy and entropy determine affinity.
**Sources**: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Basic docking — Autodock Vina 1.2.0 documentation.md; autodock-vina-readthedocs-io-en-latest.pdf; Molecular dynamics simulations and drug discovery.md; Protein–Ligand Docking in the Machine-Learning Era.md; Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md
**Last updated**: 2026-11-01

---

## The four principal interaction types

Molecular docking algorithms describe complexes held together non-covalently. In biological systems the four main classes are hydrogen bonds, ionic bonds, van der Waals interactions, and hydrophobic contacts (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Non-covalent interactions are individually weak — roughly 1–5 kcal/mol — compared with covalent bonds, but many act together at an interface to produce stable and specific association (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

### Hydrogen bonds

A hydrogen bond has the form D–H···A, where D is an electronegative donor atom, H is the hydrogen attached to it, and A is the acceptor. The interaction is primarily electrostatic and directional. Its strength is about **5 kcal/mol**, compared with about **110 kcal/mol** for an O–H covalent bond (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Because both partners are solvated before binding, making a protein–ligand hydrogen bond usually requires breaking protein–water and ligand–water hydrogen bonds. The net contribution to affinity is therefore much smaller than the raw bond energy and depends on the enthalpy and entropy of desolvation (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Hydrogen bonds also stabilize protein and nucleic acid secondary structure, so a ligand that competes with an internal hydrogen bond must pay that cost.

### Ionic interactions

Ionic interactions are electrostatic attractions between oppositely charged groups — for example, a protonated amine with an aspartate or glutamate carboxylate, or a deprotonated acid with a lysine ammonium. In water, ions are surrounded by hydration shells, so the strength of a salt bridge inside a buried site is much greater than on the solvent-exposed surface (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Charge–charge interactions are long-ranged and can accelerate association; they can also produce steep, orientation-sensitive scoring terms that are sensitive to protonation assignment.

### Van der Waals interactions

Transient dipoles in the electron clouds of nearby atoms produce nonspecific van der Waals forces. A typical strength is about **1 kcal/mol**, weaker than a hydrogen bond, but the large number of contacts in a well-packed interface makes the total contribution substantial (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Van der Waals interactions are modeled with the Lennard-Jones 6–12 potential: a steep repulsive r⁻¹² term and an attractive r⁻⁶ term (source: Molecular dynamics simulations and drug discovery.md). Shape complementarity — the quality of packing — is largely a van der Waals phenomenon.

### Hydrophobic interactions

Nonpolar groups tend to exclude water and aggregate. The classical explanation is an entropy gain from releasing ordered water at the solute surface, although the molecular mechanism remains debated; scaled-particle theory treats the effect as size-dependent (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Hydrophobic burial is often the dominant driving force for binding of drug-like ligands, and it is the reason the lock-and-key model can be entropy-dominated (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Additional interactions

- **CH–π interactions** occur between a polarized C–H bond and an aromatic ring and are common in protein–sugar complexes (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).
- **π–π stacking** between aromatic rings contributes to packing and orientation.
- **Halogen bonds** between a σ-hole on a halogen and a Lewis base are directional and can improve potency.
- **Metal coordination** — zinc, magnesium, calcium, and iron sites require explicit treatment; AutoDock4Zn adds tetrahedral zinc pseudo-atoms and a calibrated potential derived from 292 crystal complexes (source: autodock-vina-readthedocs-io-en-latest.pdf).
- **Water-mediated contacts** — bridging waters extend the interaction network and can be decisive for affinity (source: autodock-vina-readthedocs-io-en-latest.pdf; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## From interactions to affinity

The binding free energy is the balance of enthalpy and entropy:

$$\Delta G = \Delta H - T\Delta S$$

Both terms are global properties of the protein–ligand–solvent system, not sums of isolated contacts (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Enthalpy

Binding enthalpy reflects the net heat of forming and breaking interactions: protein–solvent and ligand–solvent contacts are lost, protein–ligand contacts and solvent–solvent contacts are gained, and solvent reorganizes around the complex. Each contribution can be favorable or unfavorable, and the net is their sum (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

A large negative ΔH usually means many strong polar contacts, but it is almost always accompanied by an entropic penalty from restricted motion — the basis of enthalpy–entropy compensation (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Entropy

The binding entropy has three separable contributions (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

$$\Delta S = \Delta S_{solv} + \Delta S_{conf} + \Delta S_{r/t}$$

- **ΔS_solv** — solvent entropy change from surface burial and water release; usually positive and favorable.
- **ΔS_conf** — conformational entropy change of protein and ligand; can be favorable or unfavorable.
- **ΔS_r/t** — loss of translational and rotational degrees of freedom when two particles become one; unfavorable.

Binding must overcome the unavoidable ΔS_r/t penalty through solvent entropy gain, favorable enthalpy, or both (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## How docking treats interactions

Docking scoring functions approximate these non-covalent terms with highly simplified potentials:

- **Vina** uses a united-atom empirical function combining steric attraction (Gaussian terms), quadratic repulsion, hydrophobic contacts, and directional hydrogen bonding, coupled with an explicit ligand torsional penalty term ($N_{\text{rot}}$); atomic partial charges are ignored during score evaluation (source: Basic docking — Autodock Vina 1.2.0 documentation.md; autodock-vina-readthedocs-io-en-latest.pdf).
- **AutoDock4** adds an explicit Coulomb term and uses Gasteiger–Marsili charges (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Knowledge-based functions** derive distance-dependent potentials from structural databases and therefore implicitly encode interaction preferences (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md).

The gap between these static approximations and real thermodynamics is the main reason docking scores correlate poorly with measured affinities (source: Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md). Furthermore, as demonstrated in molecular simulation benchmarks, drug unbinding kinetics ($k_{\text{off}}$) and residence times are frequently controlled by specific metastable hydrogen-bond networks or bridging waters along transition pathways rather than overall contact numbers (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[benchmark-drug-targets]].

## Interaction fingerprints

Because affinity is hard to score but interaction patterns are easier to recognize, post-docking filters convert the 3D contact pattern into a bit string — an interaction fingerprint (IFP) — and compare it with a reference co-crystal pose. IFP filtering achieved the best early enrichment on LIT-PCBA (average EF1% = 7.46), outperforming classical and ML scoring functions on that benchmark (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Related pages
- [[molecular-recognition-models]]
- [[binding-thermodynamics]]
- [[binding-kinetics]]
- [[water-and-solvation]]
- [[scoring-functions]]
- [[protein-ligand-docking]]
- [[benchmark-drug-targets]]
- [[du-2016-protein-ligand-interactions]]
- [[biochemistry-free-for-all]]
