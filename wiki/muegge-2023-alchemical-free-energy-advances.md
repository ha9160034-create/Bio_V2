# Muegge & Hu (2023) — Alchemical Free Energy in Drug Discovery

**Summary**: An industry perspective detailing recent advances in relative and absolute binding free energy calculations, covering fragment-to-lead growth, homology models, water sampling, and covalent inhibitors.
**Sources**: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Ingo Muegge and Yuan Hu, "Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery," *ACS Medicinal Chemistry Letters* 14(3):244–250 (2023), doi:10.1021/acsmedchemlett.2c00541 (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).

The review assesses the practical translation of Alchemical Binding Free Energy (BFE) methodologies into prospective industrial drug discovery pipelines. It contrasts Relative Binding Free Energy (RBFE) with Absolute Binding Free Energy (ABFE) and examines recent breakthroughs in applying alchemical methods to fragment optimization, homology models, buried water networks, and covalent inhibitors (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).

## RBFE versus ABFE in prospective workflows

The authors contrast the two primary alchemical calculation paradigms (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf):
- **Relative Binding Free Energy (RBFE)**: Mutates one ligand into another along a congeneric series using alchemical thermodynamic cycles. In prospective medicinal chemistry campaigns, RBFE routinely achieves an accuracy of **$\sim 1$ kcal/mol (RMSE)**, driving lead optimization by ranking prospective designs prior to chemical synthesis.
- **Absolute Binding Free Energy (ABFE)**: Alchemically decouples an entire ligand from the receptor pocket and from bulk solvent into non-interacting ghost states. Requiring orientational/positional Boresch restraints and standard state volume corrections ($1\text{ M} = 1,660 \text{ \AA}^3/\text{molecule}$), ABFE calculates absolute $\Delta G_{\text{bind}}$ across non-congeneric, unrelated chemotypes, albeit requiring 5- to 10-fold greater sampling effort than RBFE (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).

## Key advancements in challenging applications

### 1. Fragment-to-lead progression
Fragment-based drug design (FBDD) relies on low-affinity (micromolar to millimolar) starting points that require systematic growing, linking, or merging. In a systematic benchmark evaluating 90 fragments across eight protein systems, alchemical BFE accurately predicted favorable growth vectors and side-chain elaborations, substantially accelerating hit-to-lead timelines (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).

### 2. Homology and machine-learning models
Historically, alchemical BFE was restricted to high-resolution (<2.0 Å) crystal structures. The authors highlight a comprehensive study evaluating congeneric ligand series in 29 homology models across 14 protein targets:
- Models derived from templates sharing **30%–50% sequence identity** yielded FEP prediction accuracies that were **comparable to crystal structure baselines**, provided that the initial pocket was subjected to induced-fit MD relaxation prior to perturbation (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).
- This finding expands the domain of applicability of alchemical free energy to homology models and AlphaFold-predicted structures where experimental holo complexes are absent.

### 3. Buried waters and Grand Canonical Monte Carlo
In deep, hydrophobic cavities, trapped structural water molecules display slow exchange rates (>microseconds) with bulk solvent, leading to severe hysteresis and perturbation failures when substituents displace or trap waters. Modern AFE engines integrate Grand Canonical Monte Carlo (GCMC) or WaterMap algorithms into the alchemical simulation, dynamically inserting or deleting water molecules as the ligand changes shape to maintain true thermodynamic equilibrium (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).

### 4. Charged perturbations
Alchemical mutations that involve a net charge change (e.g., neutral ligand to anion/cation) create strong finite-size electrostatic artifacts in periodic boundary conditions under Particle Mesh Ewald (PME). Modern workflows eliminate these artifacts by applying analytical Rocklin/Hummer finite-size corrections or by employing co-alchemical counter-ion transformations in bulk solvent (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).

### 5. Covalent inhibitors
For targeted covalent therapeutics, total potency is governed by a two-step mechanism: initial non-covalent association ($K_i$) followed by chemical bond formation ($k_{\text{inact}}$), characterized by the efficiency constant $k_{\text{inact}}/K_i$. The authors emphasize that alchemical BFE calculations are used to optimize the pre-reaction non-covalent binding affinity ($K_i$), ensuring correct orientation of the reactive warhead prior to nucleophilic attack (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).

## Industrial impact and cycle times

The authors note that prospective deployment of FEP prioritizes 20 to 50 virtual designs per synthesis cycle, decreasing the number of required compound syntheses by 30%–50% while accelerating the arrival at preclinical drug candidates (source: Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf).

## Related pages
- [[binding-free-energy]]
- [[structure-based-drug-design]]
- [[molecular-dynamics-simulation]]
- [[water-and-solvation]]
- [[force-fields]]
- [[virtual-screening]]

