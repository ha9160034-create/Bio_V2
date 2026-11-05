# Poli et al. (2020) — MM-PBSA Methods in Virtual Screening

**Summary**: A review evaluating the practical protocols, enrichment improvements, and parameter calibration for applying MM-PBSA and MM-GBSA rescoring in structure-based virtual screening.
**Sources**: Application of MM-PBSA Methods in Virtual Screening.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Giulio Poli, Carlotta Granchi, Flavio Rizzolio, and Tiziano Tuccinardi, "Application of MM-PBSA Methods in Virtual Screening," *Molecules* 25(8):1971 (2020), doi:10.3390/molecules25081971 (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

The review investigates the application of Molecular Mechanics/Poisson–Boltzmann Surface Area (MM-PBSA) and Molecular Mechanics/Generalized Born Surface Area (MM-GBSA) as intermediate-tier rescoring filters to rescue true active compounds and filter out false-positive docking hits in structure-based virtual screening (SBVS) (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

## The virtual screening workflow hierarchy

The authors delineate the role of MM-PBSA within a tiered, hierarchical virtual screening cascade (source: Application of MM-PBSA Methods in Virtual Screening.pdf):
1. **Primary docking screen**: High-throughput docking (e.g., Glide, AutoDock Vina, GOLD, FRED) rapidly filters chemical libraries of $10^5$ to $10^7$ molecules down to the top 500 to 2,000 candidate poses.
2. **Intermediate MM-PBSA/GBSA filter**: The selected subset is rescored with MM-PBSA/GBSA to incorporate dynamic solvent relaxation, continuum electrostatics, and Lennard-Jones van der Waals terms.
3. **Hit selection & experimental testing**: Top-rescored molecules (typically 50–200 compounds) proceed to visual inspection, structural clustering, and biochemical assay validation (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

## Rescoring protocols and computational cost

The review evaluates three principal operational protocols for MM-PBSA rescoring (source: Application of MM-PBSA Methods in Virtual Screening.pdf):

### 1. Minimization-only protocol
The docked complex undergoes conjugate-gradient energy minimization (in vacuo or in implicit solvent) followed by static MM-PBSA evaluation. While highly economical (seconds per ligand), it omits thermal fluctuations, conformational entropy, and solvent reorganization (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

### 2. Short molecular dynamics (MD) trajectory protocol
The complex is solvated in explicit water (e.g., TIP3P) with counter-ions, equilibrated, and simulated for 0.5 to 5.0 ns. Free energies are computed over 10 to 100 snapshots extracted along the trajectory. This captures induced-fit side-chain relaxation, hydrogen-bond fluctuations, and water displacement, producing significantly higher enrichment than minimization alone (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

### 3. Single-trajectory enforcement
In virtual screening applications, the 1-trajectory approach is universally utilized: snapshots of the receptor and ligand are extracted directly from the complex trajectory, which guarantees cancellation of internal bonded terms ($\Delta E_{\text{bonded}} = 0$) and minimizes statistical variance (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

## Target benchmarks and quantitative enrichment gains

The authors review benchmark studies across diverse enzyme and receptor families (source: Application of MM-PBSA Methods in Virtual Screening.pdf):
- **COX-2, Estrogen Receptor, and p38 MAP Kinase**: Rescoring top docking hits (from FRED/ChemScore or Glide) with MD-based MM-PBSA increased 1% Enrichment Factors (EF1%) from 1.2–2.0 to >5.0 (source: Application of MM-PBSA Methods in Virtual Screening.pdf).
- **Proteases (Thrombin, BACE1, HIV-1 Protease)**: False-positive rates decreased by over 40% after MM-PBSA rescoring, primarily because poses that established superficial docking contacts dissociated or showed severe polar desolvation penalties during MD simulation (source: Application of MM-PBSA Methods in Virtual Screening.pdf).
- **Neuraminidase and Metalloenzymes**: The inclusion of explicit structural water molecules and the adjustment of solute dielectric constants ($\epsilon_{\text{int}} = 4$) were required to prevent electrostatic over-polarization in charged active sites (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

## Parameter calibration in virtual screening

The authors highlight three critical practical guidelines for VS practitioners (source: Application of MM-PBSA Methods in Virtual Screening.pdf):
1. **Solute dielectric constant ($\epsilon_{\text{int}}$)**: Standard hydrophobic cavities perform optimally with $\epsilon_{\text{int}} = 1$ or $2$. Highly charged or solvent-accessible binding clefts require $\epsilon_{\text{int}} = 4$ (or even up to 10) to avoid overestimating electrostatic interaction energies.
2. **Entropy omission in congeneric ranking**: Computing the normal mode vibrational entropy term ($-T\Delta S$) is computationally prohibitive in virtual screening and frequently degrades enrichment due to high numerical noise. Omitting $-T\Delta S$ saves >80% of calculation time while preserving or improving ranking accuracy.
3. **Structural water molecules**: Water molecules that mediate ligand–protein bridges and display long residence times in MD simulations should be treated as part of the rigid receptor during PB continuum boundary definition (source: Application of MM-PBSA Methods in Virtual Screening.pdf).

## Related pages
- [[virtual-screening]]
- [[binding-free-energy]]
- [[protein-ligand-docking]]
- [[molecular-dynamics-simulation]]
- [[scoring-functions]]
- [[water-and-solvation]]

