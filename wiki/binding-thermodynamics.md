# Binding Thermodynamics

**Summary**: Binding thermodynamics decomposes affinity into Gibbs free energy, enthalpy, and entropy, and connects them to the binding constant; isothermal titration calorimetry is the only method that measures binding enthalpy directly, while SPR and fluorescence polarization provide complementary affinity and kinetic data.
**Sources**: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md
**Last updated**: 2026-11-01

---

## The core equations

For the reversible reaction P + L ⇌ PL at equilibrium (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

$$K_b = \frac{[PL]}{[P][L]} = \frac{k_{on}}{k_{off}} = \frac{1}{K_d}$$

with K_b in M⁻¹ and K_d in M. The standard free energy is:

$$\Delta G^\circ = -RT\ln K_b$$

and the general free energy at any reaction quotient Q is:

$$\Delta G = \Delta G^\circ + RT\ln Q$$

When Q = K_b, ΔG = 0. The enthalpic and entropic contributions are:

$$\Delta G = \Delta H - T\Delta S$$

Binding occurs spontaneously only when ΔG is negative; the magnitude of the negative ΔG determines complex stability, or affinity (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Enthalpy

Enthalpy is the total energy of the system plus the PV work needed to make room for it. ΔH is negative for exothermic processes (formation of favorable non-covalent interactions) and positive for endothermic ones (disruption of favorable interactions) (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The binding enthalpy is often loosely described as the energy of the non-covalent interactions formed at the interface. In reality it is a global property of the protein–ligand–solvent system (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- loss of protein–water and ligand–water hydrogen bonds and van der Waals contacts;
- formation of protein–ligand contacts;
- solvent reorganization near the complex;
- changes in intramolecular interactions.

Each term can be favorable or unfavorable, and the net is their sum. A large negative ΔH is not necessarily achievable without a compensating entropy loss.

## Entropy

Entropy measures how heat is distributed and, loosely, the disorder or freedom of the system. The binding entropy separates into three terms (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

$$\Delta S = \Delta S_{solv} + \Delta S_{conf} + \Delta S_{r/t}$$

- **ΔS_solv** — from surface burial and release of ordered water; usually positive and favorable because the released waters gain translational freedom.
- **ΔS_conf** — change in conformational freedom of protein and ligand; can be favorable or unfavorable.
- **ΔS_r/t** — loss of translational and rotational degrees of freedom when two molecules form one complex; always unfavorable.

Binding must overcome the ΔS_r/t penalty through solvent entropy gain, favorable enthalpy, or both (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Enthalpy–entropy compensation

Strong binding from many favorable interactions produces a large negative ΔH but also restricts motion, producing a negative ΔS. Conversely, a large entropy gain usually requires breaking favorable interactions and carries an enthalpic penalty. The result is that changes in ΔH and ΔS often cancel, leaving only a modest change in ΔG — **enthalpy–entropy compensation** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The phenomenon has been debated for decades. Criticisms include (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- it may be an artifact of a narrow temperature range;
- it may arise from random or systematic experimental errors;
- it may reflect data-selection bias.

Nevertheless, it is frequently observed in calorimetric data and theoretical studies, and is generally regarded as a genuine physical phenomenon when stringent criteria are applied. It is rooted in the formation and disruption of weak non-covalent interactions and is influenced by solvent structure, pocket flexibility, ligand structure, and changes in intermolecular forces.

The practical consequence for drug design is that optimizing ΔG requires maximizing the favorable enthalpy or entropy while minimizing the compensating penalty, not simply adding interactions (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Experimental methods

### Isothermal titration calorimetry (ITC)

ITC is the only method that directly measures the heat exchanged during complex formation at constant temperature, and is the gold standard for binding thermodynamics (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

**How it works:**

1. Ligand is titrated into a solution of protein in the sample cell.
2. The instrument maintains the sample and reference cells at the same temperature by modulating power; exothermic binding decreases the power, endothermic binding increases it.
3. The primary data are power versus time, consisting of peaks whose areas are the heat per injection.
4. As the binding site saturates, peak areas decrease until only dilution heat remains.
5. The binding curve — heat per injection versus [L]/[P] — is fitted to a binding model to obtain K_b, ΔH, and stoichiometry n.
6. ΔG° is calculated from K_b, and ΔS from ΔG = ΔH − TΔS.

**What it provides:** K_b, ΔG, ΔH, ΔS, n from a single titration; ΔC_p from ΔH measured at several temperatures (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

**Advantages:** no immobilization, modification, or labelling; no molecular-weight limit; high precision (binding constant error typically ~5%); broad affinity range (~10⁻²–10³ µM, with high-affinity measurements up to 10⁸–10⁹ M⁻¹); small heat effects (0.1 µcal) and heat rates (0.1 µcal/s) detectable (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

**Disadvantages:** heat is a universal signal, so binding heat must be separated from dilution, mixing, and temperature artifacts; low signal when ΔH is small; slow processes may be missed; large sample consumption; traditionally low throughput, although robotic and array-based nanocalorimeters are changing this (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Surface plasmon resonance (SPR)

SPR measures changes in refractive index near a gold sensor surface. The protein is immobilized, ligand is injected, and the response in resonance units rises during association and falls during dissociation. Fitting the sensorgram yields k_on and k_off, and K_b = k_on/k_off (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

- **Strength:** real-time kinetics, label-free, higher affinity range than ITC (typically 10⁻⁶–10 µM) (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Limitation:** mass-transport limitation makes k_on faster than ~10⁶ M⁻¹·s⁻¹ difficult to measure accurately; immobilization perturbs conformational and translational/rotational entropy and therefore the association rate (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Enthalpy:** van't Hoff analysis of K_d versus temperature gives ΔH with much less protein than ITC, though less rigorously (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Fluorescence polarization (FP)

FP exploits the fact that polarized emission from a fluorescently labelled ligand depolarizes faster when the ligand is free than when bound. Competition experiments displace the labelled ligand with unlabelled compounds and yield IC₅₀, converted to K_i or K_d by the Cheng–Prusoff equation. van't Hoff analysis of K_d versus temperature gives ΔH (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

- **Advantages:** mix-and-read, no separation steps, low reagent consumption, inexpensive equipment, samples can be re-read, easily automated for HTS.
- **Limitations:** the signal is proportional rather than direct, so affinity can depend on assay conditions; autofluorescence, light scattering, label effects, non-specific binding, and aggregation can confound results.

### Structural and dynamic methods

X-ray crystallography, NMR, Laue diffraction, SAXS, cryo-EM, single-molecule fluorescence, and time-resolved hydrogen–deuterium exchange MS provide structural and dynamic information that complements thermodynamics (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). NMR is particularly valuable because it covers picoseconds to seconds and can probe entropic contributions.

## Heat capacity

The change in heat capacity ΔC_p is obtained from the slope of ΔH versus temperature:

$$\Delta H(T) = \Delta H(T_0) + \Delta C_p(T-T_0)$$

ΔC_p correlates strongly with buried surface area. Burial of apolar surface gives negative ΔC_p; burial of polar surface gives positive ΔC_p. Because the hydration of hydrophobic surfaces has a different heat capacity from bulk water, ΔC_p links thermodynamic parameters to structural information (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Worked example: HIV-1 protease inhibitor optimization

Freire and coworkers used ITC to optimize HIV-1 protease inhibitors, comparing pairs that differ by a single functional group (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- **Thioether → sulfonyl (KNI-10033 → KNI-10075):** ΔH became more negative by 3.9 kcal/mol, but the entropy gain fell by 4.2 kcal/mol, so ΔG rose by 0.3 kcal/mol; a strong new hydrogen bond did not improve affinity because the entropic penalty was larger (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Methyl → hydroxyl (KNI-10052 → KNI-10054):** ΔH became more negative by 4.4 kcal/mol, overcoming a 3.9 kcal/mol entropy loss and lowering ΔG by 0.5 kcal/mol — a factor of 2 in affinity (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Hydrogen → methyl (KNI-10046 → KNI-10030):** both enthalpy (−0.8 kcal/mol) and entropy (+0.5 kcal/mol) improved, lowering ΔG from −9.6 to −10.9 kcal/mol, an 8.7-fold affinity increase, though increasing hydrophobicity (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The lesson is that a single polar group can improve enthalpy without improving affinity, and that enthalpically optimized compounds are becoming a design trend despite being harder to achieve (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Thermodynamics of the recognition models

- **Lock-and-key** — entropy-dominated; solvent entropy gain from water displacement drives binding (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Induced fit** — enthalpy-dominated; new interactions must pay for desolvation, disruption of existing interactions, and entropy loss (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Conformational selection** — entropy and enthalpy act sequentially: solvent entropy gain in selection, enthalpy decrease in adjustment (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

See [[molecular-recognition-models]].

## Computational thermodynamics

MD and free-energy methods aim to reproduce these thermodynamic quantities (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). Classical docking scoring functions attempt to approximate the master free energy equation ($\Delta G = \Delta H - T\Delta S$) through empirical, knowledge-based, or physics-based terms representing electrostatics, van der Waals interactions, hydrogen bonding, desolvation, and ligand torsional entropy losses (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). However, rigid docking intrinsically neglects protein conformational reorganization enthalpy and oversimplifies solvent entropy (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Alchemical calculations give ΔG with ~1–2 kcal/mol accuracy; the decomposition into ΔH and ΔS is much less reliable because free energy is variational and entropy is not (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). In a host–guest study, ΔG from two force fields and two routes agreed within 2 kcal/mol, but the ΔH/ΔS split depended strongly on the force field (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). High-precision experimental validation is provided by [[experimental-binding-assays]] and benchmark targets such as HIV-1 protease (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). See [[benchmark-drug-targets]].

## Related pages
- [[binding-free-energy]]
- [[binding-kinetics]]
- [[experimental-binding-assays]]
- [[protein-ligand-interactions]]
- [[molecular-recognition-models]]
- [[water-and-solvation]]
- [[force-fields]]
- [[benchmark-drug-targets]]
- [[du-2016-protein-ligand-interactions]]
