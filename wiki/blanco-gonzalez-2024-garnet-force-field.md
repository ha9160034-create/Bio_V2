# Blanco-González et al. (2024) — Garnet Force Field

**Summary**: Research from the MRC Laboratory of Molecular Biology presenting Garnet, an end-to-end graph neural network that parameterizes classical molecular mechanics force fields from scratch for both proteins and small molecules using continuous atom typing and a double exponential nonbonded potential.
**Sources**: `Training a force field for proteins and small molecules from.pdf`
**Last updated**: 2026-11-02

---

## Machine-learning force field parameterization from scratch

Blanco-González et al. introduce **Garnet**, an automated force field framework developed at the MRC Laboratory of Molecular Biology that parameterizes classical molecular dynamics (MD) simulations without relying on human-curated atom types or legacy force field parameters (source: Training a force field for proteins and small molecules from.pdf). 

For decades, biomolecular force field development has been bifurcated: protein force fields (e.g., AMBER, CHARMM) were fitted to small peptide model systems, while small-molecule general force fields (e.g., GAFF, CGenFF) were constructed via separate manual atom-typing rules (source: Training a force field for proteins and small molecules from.pdf). Even modern graph-based approaches such as Espaloma parameterize valence terms via neural networks but continue to inherit legacy van der Waals parameters (Lennard-Jones $\sigma$ and $\varepsilon$) and predefined water models (source: Training a force field for proteins and small molecules from.pdf). Garnet resolves this fragmentation by learning all bonded and nonbonded parameters—for proteins, small molecules, and water simultaneously—from quantum mechanical (QM) data and condensed-phase experimental observables (source: Training a force field for proteins and small molecules from.pdf).

## Continuous atom typing via Graph Neural Networks (GNN)

Instead of discretizing chemical space into rigid alphanumeric atom types, Garnet implements continuous atom typing:
- A molecular topology graph is processed by a message-passing **Graph Neural Network (GNN)** with a 2-bond receptive field, generating continuous 64-dimensional atom embeddings (source: Training a force field for proteins and small molecules from.pdf).
- Atom embeddings are projected via lightweight multilayer perceptrons (MLPs) to predict equilibrium bond lengths ($r_0$), bond force constants ($k_b$), equilibrium angles ($\theta_0$), angle force constants ($k_\theta$), and periodic dihedral Fourier coefficients ($k_n, \delta_n$) (source: Training a force field for proteins and small molecules from.pdf).
- Atomic partial charges are derived via a charge-equilibration layer trained to reproduce Minimal Basis Iterative Stockholder (**MBIS**) partial charges extracted from the large-scale **SPICE** dataset (source: Training a force field for proteins and small molecules from.pdf).

## The double exponential nonbonded potential

A major technical breakthrough in Garnet is the replacement of the standard 12-6 Lennard-Jones (LJ) potential with the **double exponential potential**:
- During end-to-end backpropagation through MD trajectories via ensemble gradient techniques, the steep $r^{-12}$ repulsive term of the Lennard-Jones potential caused severe numerical gradient instability and exploding optimization trajectories (source: Training a force field for proteins and small molecules from.pdf).
- Garnet adopts a double exponential functional form:
  $$U_{\text{dexp}}(r) = \varepsilon \left[ \frac{\beta}{\alpha - \beta} \exp\left(\alpha \left(1 - \frac{r}{\sigma}\right)\right) - \frac{\alpha}{\alpha - \beta} \exp\left(\beta \left(1 - \frac{r}{\sigma}\right)\right) \right]$$
  where $\sigma$ and $\varepsilon$ represent atomic size and dispersion well-depth, and $\alpha$ and $\beta$ are global parameters governing repulsive steepness and attractive decay (source: Training a force field for proteins and small molecules from.pdf).
- The exponential repulsion provides a more physically realistic representation of quantum mechanical Pauli exchange repulsion than the empirical $r^{-12}$ term, yields numerically stable derivatives during end-to-end learning, and executes with negligible GPU overhead in OpenMM via custom force expressions (source: Training a force field for proteins and small molecules from.pdf).
- The water model is co-trained alongside the solute parameters, eliminating historical mismatches between solvent models (e.g., TIP3P/TIP4P) and solute van der Waals radii (source: Training a force field for proteins and small molecules from.pdf).

## Benchmark performance across scales

Garnet was evaluated across a diverse hierarchy of molecular benchmarks:

| Benchmark Domain | Metric / Dataset | Garnet Performance | Comparison Baselines |
| :--- | :--- | :--- | :--- |
| **Quantum Mechanical Forces** | SPICE Test Set 1 v1.0 (mean force error) | **$398\text{ kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$** | Espaloma 0.4.0: $522$; OpenFF 2.2.1: $1104$; MACE-OFF23 (MLIP): $40$ |
| **Potential Energy Differences** | SPICE Test Set 1 v1.0 ($\Delta E$ error) | **$14.2\text{ kJ/mol}$** | Espaloma 0.4.0: $21.3$; OpenFF 2.2.1: $36.5$; MACE-OFF23: $3.3$ |
| **Folded Proteins** | 1-µs simulations of globular proteins (GB3, Ubiquitin) | Stable trajectories; low backbone RMSD ($< 1.5\text{ \AA}$) | Comparable to Amber14SB/TIP3P |
| **Intrinsically Disordered Proteins (IDPs)** | Radius of gyration ($R_g$) and NMR chemical shifts | Captures expanded conformational ensembles without over-compaction | Outperforms classical force fields prone to hydrophobic collapse |
| **Protein–Ligand Relative Binding Free Energy (RBFE)** | Benchmark targets (CDK2, Tyk2, PTP1B, MCL1) | High correlation with experiment ($\text{RMSE} \approx 1.1\text{--}1.4\text{ kcal/mol}$) | Competitive with Amber14SB + OpenFF 2.2.0 and OPLS4 |

(source: Training a force field for proteins and small molecules from.pdf)

## Implications for computer-aided drug discovery

Garnet demonstrates that an end-to-end machine-learning architecture can discover a unified, fully differentiable force field for both proteins and small molecules from raw physical data (source: Training a force field for proteins and small molecules from.pdf). By circumventing discrete atom typing and manual parameter tuning, the framework enables rapid continuous adaptation to non-standard ligands, covalent modifications, non-canonical amino acids, and metalloenzyme environments in structure-based drug design (source: Training a force field for proteins and small molecules from.pdf).

## Related pages
- [[open-force-field-and-ligand-parameterization]]
- [[wang-2024-open-force-field-initiative]]
- [[force-fields]]
- [[binding-free-energy]]
- [[alchemical-free-energy-methods]]
- [[simulation-and-benchmark-best-practices]]
- [[ligand-preparation-and-protonation]]
- [[molecular-dynamics-simulation]]

