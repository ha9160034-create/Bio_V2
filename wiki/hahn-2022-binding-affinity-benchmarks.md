# Hahn et al. (2022) — Protein–Ligand Binding Affinity Benchmarks

**Summary**: The LiveCoMS comprehensive best practices guide for constructing, preparing, and evaluating protein–ligand binding affinity benchmarks, introducing the PLBenchmarks repository and the arsenic statistical assessment toolkit.
**Sources**: `Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf`
**Last updated**: 2026-11-02

---

## Benchmark construction and the domain of applicability

Hahn et al. provide authoritative standards for establishing meaningful benchmarks in computer-aided drug discovery (CADD), specifically targeting relative binding free energy (RBFE) and absolute binding free energy (ABFE) methods (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). The authors observe that retrospective free-energy literature often suffers from poorly curated experimental data, ambiguous protonation states, uncontrolled structural artifacts, and statistically deficient validation, leading to benchmark results that fail to predict prospective performance in lead optimization (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

To address these limitations, the authors introduce **PLBenchmarks**, an open-access, versioned, and standardized repository of curated protein–ligand benchmark datasets, alongside **arsenic**, an open-source Python toolkit implementing standardized best practices for evaluating free-energy calculations (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

## Curation of experimental bioactivity data

The quality of any computational benchmark is strictly bounded by the fidelity of its experimental reference data:
- **Assay Type Hierarchies**: Direct thermodynamic and biophysical binding assays—specifically Isothermal Titration Calorimetry (ITC) and Surface Plasmon Resonance (SPR)—are the gold standard because they directly monitor equilibrium binding without enzymatic coupling (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Well-controlled fluorescence polarization (FP) or competitive radioligand displacement assays with verified equilibrium conditions represent acceptable secondary sources (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **The Problem with $\text{IC}_{50}$ Values**: Functional inhibition values ($\text{IC}_{50}$) derived from cell-based or enzymatic activity assays should be avoided where possible (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Converting $\text{IC}_{50}$ to inhibition constants ($K_i$) via the classical Cheng–Prusoff relationship ($K_i = \text{IC}_{50} / (1 + [S]/K_m)$) requires accurate knowledge of substrate concentration $[S]$ and the Michaelis constant $K_m$; variations in assay conditions, enzyme pre-incubation, or tight-binding regimes ($[E] \approx \text{IC}_{50}$) introduce non-thermodynamic shifts that invalidate free-energy comparisons (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Experimental Dynamic Range**: Benchmarks require an experimental affinity dynamic range spanning at least $2\text{ to }3$ orders of magnitude (or $\Delta\Delta G \ge 2.7\text{--}4.1\text{ kcal/mol}$) to provide sufficient statistical power for correlation analysis (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Experimental Uncertainty Floor**: Replicate biophysical measurements typically exhibit an inherent experimental noise floor of $\sim 0.3\text{ to }0.5\text{ kcal/mol}$ (or roughly a factor of $2$ in $K_d$); computational models achieving root-mean-square errors (RMSE) within this window cannot be further differentiated without higher-precision experimental techniques (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

## Structural quality criteria and system preparation

Rigorous structural vetting of protein–ligand complexes is mandatory:
- **Resolution and Refinement Metrics**: High-resolution X-ray crystal structures ($\le 2.0\text{--}2.5\text{ \AA}$) with low free $R$-factors ($R_{\text{free}} < 0.25$) and high Real-Space Correlation Coefficients (RSCC) are preferred (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Electron density difference maps ($2F_o - F_c$ and $F_o - F_c$) must be visually inspected to confirm unambiguous placement of the ligand scaffold, substituents, and binding-site sidechains (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Multiple Crystal Forms and Packing Artifacts**: Crystal lattice contacts adjacent to the binding pocket must be checked to ensure that crystal packing forces have not distorted active-site loops into unphysiological conformations (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Protonation and Tautomeric States**: The protonation states of ionizable amino acids (His, Asp, Glu, Lys, Cys) and ligand tautomers must be evaluated at the specific pH of the bioactivity assay (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Shifts in $\text{p}K_a$ upon binding must be accounted for when relating thermodynamic cycles to experimental observables (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Cofactors and Structural Water**: Crystallographic waters that bridge ligand–protein hydrogen-bond networks or mediate catalytic metal coordination (e.g., zinc or magnesium) must be explicitly retained, while disordered bulk solvent is removed (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

## Perturbation network design for relative free energy (RBFE)

When benchmarking alchemical transformations across a congeneric ligand series:
- **Perturbation Graph Topology**: Redundant, closed-loop perturbation networks are mandatory (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Closed thermodynamic cycles allow calculation of cycle closure hysteresis ($\sum_{\text{cycle}} \Delta\Delta G = 0$), serving as an intrinsic diagnostic of simulation convergence and phase-space sampling adequacy (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Maximum Common Substructure (MCS)**: Ligand pairs should share a large common core to minimize the number of perturbed non-hydrogen atoms (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf). Transformations involving major scaffold hops, core ring closures, or substantial charge changes introduce phase-space overlap bottlenecks that require specialized alchemical staging or absolute binding free energy (ABFE) workflows (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

## Statistical evaluation best practices

Hahn et al. outline rigorous statistical practices for benchmarking:
- **Avoid Over-Reliance on $R^2$**: The coefficient of determination ($R^2$) is highly sensitive to the experimental dynamic range; a dataset with artificially wide affinity spans can yield high $R^2$ despite severe localized prediction errors (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Recommended Metrics**: Studies should report root-mean-square error (RMSE), mean absolute error (MAE / MUE), and rank-order correlations such as Spearman's $\rho$ and Kendall's $\tau$ (which evaluates pairwise concordant vs. discordant ranking) (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).
- **Confidence Intervals**: All metrics must be accompanied by non-parametric bootstrapped $95\%$ confidence intervals to prevent overinterpreting marginal performance differences between force fields or sampling algorithms (source: Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf).

## Related pages
- [[simulation-and-benchmark-best-practices]]
- [[benchmarking-datasets]]
- [[binding-free-energy]]
- [[alchemical-free-energy-methods]]
- [[braun-2019-foundations-molecular-simulation]]
- [[open-force-field-and-ligand-parameterization]]
- [[receptor-preparation]]
- [[ligand-preparation-and-protonation]]

