# Du et al. (2016) — Insights into Protein–Ligand Interactions

**Summary**: A review that unifies the thermodynamics, kinetics, and molecular-recognition models of protein–ligand binding with the experimental methods (ITC, SPR, fluorescence polarization) and computational methods (docking, free-energy calculations) used to measure and predict them.
**Sources**: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md
**Last updated**: 2026-11-01

---

## Citation and scope

Xing Du, Yi Li, Yuan-Ling Xia, Shi-Meng Ai, Jing Liang, Peng Sang, Xing-Lai Ji, and Shu-Qun Liu, "Insights into Protein–Ligand Interactions: Mechanisms, Models, and Methods," *International Journal of Molecular Sciences* 17(2):144 (2016), doi:10.3390/ijms17020144. PMC4783901 (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The review is unusual in joining three literatures that are often treated separately: the physical chemistry of binding (kinetics and thermodynamics), the conceptual models of molecular recognition, and the experimental and computational methods used to study them. It is the wiki's main source for [[binding-thermodynamics]], [[binding-kinetics]], [[molecular-recognition-models]], and [[protein-ligand-interactions]].

## Structure of the review

1. Binding kinetics.
2. Basic thermodynamic concepts and relationships.
3. Binding driving forces and enthalpy–entropy compensation.
4. Binding models: lock-and-key, induced fit, conformational selection, and their relationships.
5. Experimental methods: ITC, SPR, fluorescence polarization, and structural methods.
6. Computational methods: docking and free-energy calculations.
7. Conclusions.

## Binding kinetics

For P + L ⇌ PL, the time-dependent association is governed by k_on and k_off. At equilibrium K_b = k_on/k_off = 1/K_d, so the kinetic constants determine the thermodynamic stability of the complex (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). See [[binding-kinetics]].

## Thermodynamics

The standard free energy is ΔG° = −RT ln K_b, with R = 1.987 cal·K⁻¹·mol⁻¹, at 1 atm, 298 K, and 1 M effective reactant concentrations. At any moment ΔG = ΔG° + RT ln Q, and at equilibrium Q = K_b so ΔG = 0. The enthalpic and entropic contributions satisfy ΔG = ΔH − TΔS (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The binding entropy separates into three terms:

$$\Delta S = \Delta S_{solv} + \Delta S_{conf} + \Delta S_{r/t}$$

- **ΔS_solv** — solvent entropy gain from surface burial and water release; usually positive and favorable.
- **ΔS_conf** — change in conformational freedom of protein and ligand; may be favorable or unfavorable.
- **ΔS_r/t** — loss of translational and rotational degrees of freedom when two particles form one complex; always unfavorable.

Binding must overcome the inevitable ΔS_r/t penalty through a large solvent entropy gain, favorable enthalpy, or both (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Enthalpy is global, not local

The binding enthalpy is often described as the energy of the non-covalent interactions formed at the interface. In reality it is a global property of the protein–ligand–solvent system: it includes loss of protein–water and ligand–water hydrogen bonds and van der Waals contacts, formation of protein–ligand contacts, solvent reorganization near the complex, and changes in intramolecular interactions. Each contribution may be favorable or unfavorable (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Enthalpy–entropy compensation

Strong binding from many favorable interactions lowers ΔH but restricts motion, lowering ΔS; a large entropy gain usually requires breaking favorable interactions and carries an enthalpic penalty. The result is a modest change in ΔG — **enthalpy–entropy compensation** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Criticisms are that the effect can be an artifact of a narrow temperature or free-energy range, of random or systematic experimental error, or of data-selection bias. Nevertheless it is frequently observed in calorimetric data and supported by theory, and is regarded as genuine when stringent criteria are applied. It is rooted in the formation and disruption of weak non-covalent interactions and is influenced by solvent structure, pocket flexibility, ligand structure, and changes in intermolecular forces (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The design implication is that optimizing affinity means maximizing the favorable enthalpy or entropy while minimizing the compensating penalty, not simply adding interactions (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). See [[binding-thermodynamics]].

## Binding models

### Diffusion and the encounter complex

Binding requires collisions that form an encounter complex. Molecular diffusion is entropy-driven and has two origins: the kinetic energy of the solute and collisions with water molecules. Long-range electrostatic attraction between oppositely charged partners can accelerate association beyond the diffusion limit. Brownian dynamics can predict diffusion-controlled association rates and encounter pathways (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Lock-and-key

Fisher's model requires rigid, perfectly matched interfaces. Studies of complementary colloidal particles show that lock-and-key association can be entropically driven by depletion forces even without attractive interactions, depending only on size, shape, and surface roughness. When the matched interfaces collide, hydration water is displaced, increasing the accessible volume and solvent entropy. Because neither partner changes conformation, ΔS_conf = 0, so the process is **entropy-dominated** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Induced fit

Koshland's model assumes a flexible binding site that rearranges around the ligand. The initial contacts must be strong enough to keep the encounter complex alive long enough for rearrangement, and some pre-existing complementarity is required. The subsequent adjustment forms new interactions at the cost of disrupting old ones, desolvation, and loss of conformational and translational/rotational entropy. The net entropy change is usually negative but is overcompensated by a larger negative enthalpy, so induced fit is **enthalpy-dominated** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Experimental examples (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- Designed receptor diastereomer 4 binding NMe4I: ΔG° = −9.1 kcal/mol, TΔS° = −0.2 kcal/mol, ΔH° = −9.3 kcal/mol — strongly enthalpy-driven.
- Barbital and phenobarbital binding to a flexible macrocyclic barbiturate receptor: six strong hydrogen bonds induce a complementary site, and the negative enthalpy overcompensates the loss of conformational entropy.
- Nitrophenyl-α-galactoside binding to the MelB permease of *Salmonella typhimurium*: ΔH = −10.3 kcal/mol, −TΔS = −3.8 kcal/mol — solely enthalpy-driven.

### Conformational selection

Conformational selection derives from free-energy-landscape theory. The protein populates an ensemble of states, and the ligand selectively binds the state whose site best matches it, shifting the population toward that state. The process has three stages (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

1. **Selective binding** — the ligand chooses a pre-existing competent conformer; resembles lock-and-key and is dominated by solvent entropy gain.
2. **Conformational adjustment** — side chains or backbone rearrange to optimize interactions; resembles induced fit and is dominated by enthalpy decrease.
3. **Population shift** — complex formation changes the barrier heights between the bound state and adjacent states, redistributing the ensemble.

In this model entropy and enthalpy act **sequentially**, not one throughout. Conformational selection also allows one protein to bind structurally distinct partners and to accommodate mutations through shifts in the landscape, which is evolutionarily advantageous (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### How the models relate

Under the funnel-like free-energy landscape (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- **Lock-and-key** is an extreme of conformational selection: a rigid protein has a smooth funnel with a single well and one conformation.
- **Conformational selection** applies to a flexible protein with multiple wells; the ligand selects one well and shifts the population.
- **Induced fit** is a key step within conformational selection, because the enhanced interactions formed after the encounter accelerate the population shift.
- Conversely, induced fit contains an implicit selection step, because the initial contacts must persist for rearrangement.

The balance depends on **ligand concentration** and **protein dynamics timescale**: conformational selection dominates at low ligand concentration and slow transitions; induced fit dominates at high concentration or fast transitions. A flux-based analysis of NADPH binding to dihydrofolate reductase and of flavodoxin folding coupled to FMN binding showed the mechanism switching from conformational selection at low ligand concentration to induced fit at high concentration (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

All three mechanisms are observed experimentally, and real binding can combine them simultaneously or sequentially (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). See [[molecular-recognition-models]].

## Experimental methods

### Isothermal titration calorimetry (ITC)

ITC is the only method that directly measures the heat of complex formation and is the gold standard for binding thermodynamics. Ligand is titrated into protein; the instrument maintains equal temperature between sample and reference cells by modulating power; the areas of the power-versus-time peaks give the heat per injection; the binding curve is fitted for K_b, ΔH, and stoichiometry n; ΔG° follows from K_b and ΔS from ΔG = ΔH − TΔS (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Advantages: no immobilization, modification, or labelling; no molecular-weight limit; high precision (K_b error typically ~5%); broad affinity range (~10⁻²–10³ µM, with high-affinity measurements up to 10⁸–10⁹ M⁻¹); detection of 0.1 µcal heat effects. Disadvantages: heat is a universal signal, so dilution and mixing heats must be separated; weak signal when ΔH is small; slow processes may be missed; large sample consumption; traditionally low throughput (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Surface plasmon resonance (SPR)

SPR measures refractive-index changes near a gold sensor surface. Protein is immobilized, ligand is injected, and the response rises during association and falls during dissociation. Fitting the sensorgram gives k_on and k_off, hence K_b. It is label-free and provides real-time kinetics, with a typical affinity range of 10⁻⁶–10 µM. Mass transport limits accurate k_on above ~10⁶ M⁻¹·s⁻¹, and immobilization perturbs the protein's conformational and translational/rotational entropy, affecting the association rate (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Fluorescence polarization (FP)

FP uses the depolarization of polarized emission from a fluorescently labelled ligand. Competition with unlabelled compounds gives IC50, converted to K_i or K_d by the Cheng–Prusoff equation. It is mix-and-read, needs no separation, consumes little reagent, and is easily automated, but the signal is proportional rather than direct, and autofluorescence, light scattering, label effects, nonspecific binding, and aggregation can confound results (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Structural and dynamic methods

X-ray crystallography, NMR, Laue diffraction, SAXS, cryo-EM, single-molecule fluorescence, and time-resolved hydrogen–deuterium exchange MS provide structural and dynamic information. NMR is particularly valuable because it covers picoseconds to seconds and can probe entropic contributions; crystallographic B-factors report some spatial distribution around the average structure (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Heat capacity

ΔC_p is obtained from the slope of ΔH versus temperature. It correlates strongly with buried surface area: apolar burial gives negative ΔC_p, polar burial gives positive ΔC_p. Because hydration water has different heat capacity from bulk water, ΔC_p links thermodynamics to structure and reports on desolvation (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Worked example: HIV-1 protease inhibitors

Freire and coworkers used ITC to compare inhibitor pairs differing by a single group (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- **Thioether → sulfonyl (KNI-10033 → KNI-10075):** ΔH became more negative by 3.9 kcal/mol, but the entropy gain fell by 4.2 kcal/mol, so ΔG rose by 0.3 kcal/mol. A strong new hydrogen bond did not improve affinity because the entropic cost was larger.
- **Methyl → hydroxyl (KNI-10052 → KNI-10054):** ΔH became more negative by 4.4 kcal/mol, overcoming a 3.9 kcal/mol entropy loss and lowering ΔG by 0.5 kcal/mol — a factor of 2 in affinity.
- **Hydrogen → methyl (KNI-10046 → KNI-10030):** both enthalpy (−0.8 kcal/mol) and entropy (+0.5 kcal/mol) improved, lowering ΔG from −9.6 to −10.9 kcal/mol, an 8.7-fold affinity increase, but the compound became more hydrophobic, with potential solubility, bioavailability, and selectivity problems.

The lesson is that a single polar group can improve enthalpy without improving affinity, and that enthalpically optimized compounds are a design trend despite being harder to achieve (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Computational methods

### Docking

Docking has two essential components: a search algorithm and a scoring function. The search space in principle includes all conformations of the unbound partners, all orientations and conformations of the ligand in a given protein state, and all protein conformations paired with all ligand conformations — far too large to explore exhaustively (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Search algorithms have evolved from rigid-body methods, through flexible-ligand methods, to flexible-ligand/flexible-protein methods. Protein flexibility is the biggest challenge: neglecting it gives poor poses, while sampling it with MD is expensive. Strategies include docking to multiple static structures (which parallels conformational selection), energy minimization with Monte Carlo or gradient descent (which parallels induced fit), rotamer search of side chains around the site, and soft-core potentials that allow small overlaps (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Scoring functions are classed as force-field-based, empirical, and knowledge-based, with consensus scoring combining several. Each family has characteristic weaknesses: force-field methods are computationally intensive and struggle with solvent and entropy; empirical methods risk over-fitting and double counting and depend on their training set; knowledge-based methods suffer from the reference-state and sparse-data problems (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). See [[protein-ligand-docking]] and [[scoring-functions]].

### Free-energy calculations

Free-energy calculations compute binding free energies from statistical thermodynamics and include energetic and entropic contributions without case-by-case parameter fitting, at a computational cost orders of magnitude above scoring functions (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The review describes three families:

- **Alchemical calculations** — the ligand is transmuted into another species or into a non-interacting dummy through nonphysical intermediates. The thermodynamic cycle gives ΔΔG_bind from the difference between the transformation in the binding site and in solution. Implementations include FEP, TI, and BAR.
- **Path sampling** — methods that compute the free-energy profile along a physical pathway.
- **Endpoint methods** — MM-PBSA/MM-GBSA and related approaches that estimate the free energy from the bound and unbound end states.

See [[binding-free-energy]] for the modern treatment and accuracy figures.

## Key numbers

| Quantity | Value |
|---|---|
| R | 1.987 cal·K⁻¹·mol⁻¹ |
| Standard state | 1 atm, 298 K, 1 M |
| ITC K_b precision | ~5% |
| ITC affinity range | ~10⁻²–10³ µM (up to 10⁸–10⁹ M⁻¹) |
| SPR k_on limit | ~10⁶ M⁻¹·s⁻¹ |
| Diastereomer 4–NMe4I | ΔG° −9.1, ΔH° −9.3, TΔS° −0.2 kcal/mol |
| MelB–nitrophenyl-α-galactoside | ΔH −10.3, −TΔS −3.8 kcal/mol |
| HIV protease KNI-10046 → KNI-10030 | ΔG −9.6 → −10.9 kcal/mol |

(source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md)

## Strengths

- Integrates thermodynamics, kinetics, and recognition models in one coherent framework.
- Gives explicit thermodynamic signatures for each recognition model, which is rare.
- Covers experimental methods in enough practical detail to be useful at the bench.
- Uses real ITC examples, especially the HIV protease series, to show why enthalpy optimization can fail.

## Limitations

- The computational sections are introductory; modern enhanced sampling, ML scoring, and ultra-large screening are not covered.
- The thermodynamic interpretation of the recognition models is presented as a rationalization, not as a proven universal rule.
- Some mechanistic conclusions (for example, the entropy dominance of lock-and-key) rely on colloidal-model experiments and analogy rather than direct protein measurements.
- Published in 2016; ITC instrumentation, SPR assays, and free-energy methods have advanced since.

## Why it matters for this wiki

This is the source for [[binding-thermodynamics]], [[binding-kinetics]], [[molecular-recognition-models]], and [[water-and-solvation]], and a supporting source for [[protein-ligand-interactions]], [[protein-ligand-docking]], and [[binding-free-energy]]. Its most distinctive contribution is the thermodynamic signature assigned to each recognition model: lock-and-key is entropy-dominated, induced fit is enthalpy-dominated, and conformational selection uses entropy and enthalpy sequentially.

## Related pages
- [[binding-thermodynamics]]
- [[binding-kinetics]]
- [[molecular-recognition-models]]
- [[protein-ligand-interactions]]
- [[water-and-solvation]]
- [[binding-free-energy]]
- [[protein-ligand-docking]]
- [[decherchi-2020-binding-thermodynamics-kinetics]]
