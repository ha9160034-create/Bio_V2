# Molecular Recognition Models

**Summary**: Three conceptual models describe how a protein recognizes a ligand — lock-and-key (rigid complementarity), induced fit (ligand-induced conformational change), and conformational selection (ligand selects a pre-existing state from a dynamic ensemble); real binding often combines them.
**Sources**: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Molecular dynamics simulations and drug discovery.md; Protein–Ligand Docking in the Machine-Learning Era.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md
**Last updated**: 2026-11-01

---

## The three models

### Lock-and-key

Fisher's lock-and-key model requires both protein and ligand to be rigid, with perfectly matched binding interfaces: only the correctly sized key enters the keyhole (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The model explains shape complementarity but cannot explain binding when the initial shapes do not match. Its thermodynamic signature is distinctive: because neither partner changes conformation, ΔS_conf is zero, and binding is driven by solvent entropy gain from displacing water at the matched interfaces. The lock-and-key process is therefore described as **entropy-dominated** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The model remains a useful first approximation for rigid pockets, and it underpins rigid-receptor docking (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### Induced fit

Koshland's induced-fit model assumes a flexible binding site that changes conformation in response to the incoming ligand, optimizing the amino-acid configuration around it — a "hand in glove" relationship (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

Induced fit is appropriate for proteins that undergo only minor conformational change on binding. Its thermodynamic signature is typically **enthalpy-dominated**: new protein–ligand interactions must be strong enough to pay for disrupting existing interactions, desolvation, and the loss of conformational and translational/rotational entropy. In the experimental examples, the net entropy change is negative but is overcompensated by a larger negative enthalpy (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Evidence includes designed host–guest systems and sugar-binding proteins:

- Receptor diastereomer 4 binding NMe4I: ΔG° = −9.1 kcal/mol, TΔS° = −0.2 kcal/mol, ΔH° = −9.3 kcal/mol (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- Nitrophenyl-α-galactoside binding to MelB permease: ΔH = −10.3 kcal/mol, −TΔS = −3.8 kcal/mol — solely enthalpy-driven (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

Induced fit requires a certain degree of pre-existing complementarity; otherwise the encounter complex is too short-lived for the rearrangement to occur (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### Conformational selection

The conformational selection model derives from free-energy-landscape theory. The native protein is not a single structure but an ensemble of conformational states with different populations. The ligand binds selectively to the state whose binding site best matches it, shifting the equilibrium toward that state (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The process has three stages (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

1. **Selective binding** — the ligand chooses a pre-existing competent conformer; this resembles lock-and-key and is dominated by solvent entropy gain.
2. **Conformational adjustment** — side chains or backbone rearrange to optimize interactions; this resembles induced fit and is dominated by enthalpy decrease.
3. **Population shift** — complex formation changes the barrier heights between the bound state and adjacent states, redistributing the ensemble toward the bound conformer.

In conformational selection, entropy and enthalpy act **sequentially** rather than one dominating throughout (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

The extended conformational selection model adds an induced-fit step after the initial selection, so the two mechanisms can coexist in one binding event (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

## How the models relate

Under the free-energy-landscape picture (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- **Lock-and-key** is an extreme of conformational selection: a very rigid protein has a smooth funnel with a single free-energy well and one conformation.
- **Conformational selection** applies to a flexible protein with multiple wells; the ligand selects one well and then shifts the population.
- **Induced fit** is a key step within conformational selection: the enhanced interactions formed after the initial encounter accelerate the population shift, extending and optimizing selection.
- Conversely, induced fit contains an implicit selection step, because the initial contacts must persist long enough for rearrangement.

The balance between conformational selection and induced fit depends on **ligand concentration** and **protein dynamics timescale** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):

- conformational selection dominates at low ligand concentration and slow conformational transitions;
- induced fit dominates at high ligand concentration or fast transitions.

A flux-based analysis of NADPH binding to dihydrofolate reductase and of flavodoxin folding coupled to FMN binding showed the mechanism switching from conformational selection at low ligand concentration to induced fit at high concentration (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

All three mechanisms are observed experimentally, and real binding can involve them simultaneously or sequentially across a broad spectrum (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Diffusion and the encounter complex

Before any recognition model applies, protein and ligand must collide. Molecular diffusion — driven by thermal energy and by collisions with water — is entropy-driven. Long-range electrostatic attraction between oppositely charged partners can accelerate association beyond the diffusion limit. Brownian dynamics can predict diffusion-controlled association rates and encounter pathways (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Implications for docking and simulation

- **Rigid-receptor docking** implicitly assumes lock-and-key. It works best for rigid pockets and fails when the competent conformation is a low-population state (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Ensemble docking** and the **relaxed complex scheme** approximate conformational selection by docking against multiple receptor conformations (source: Molecular dynamics simulations and drug discovery.md).
- **Induced fit docking** and **soft docking** approximate the induced-fit branch by allowing side-chain rearrangement or partial clashes (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Long MD and Markov state models** can distinguish the pathways directly; simulations of peptide binding to the PDZ domain and of choline binding to ChoX found both mechanisms, with real cases described as a superposition of the two idealized extremes (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- **Cryptic pockets** are a practical consequence of conformational selection: the druggable state exists in the ensemble but not in the crystal structure (source: Molecular dynamics simulations and drug discovery.md).

## Related pages
- [[protein-ligand-interactions]]
- [[protein-ligand-docking]]
- [[molecular-dynamics-simulation]]
- [[cryptic-and-allosteric-pockets]]
- [[binding-thermodynamics]]
- [[binding-kinetics]]
- [[du-2016-protein-ligand-interactions]]
