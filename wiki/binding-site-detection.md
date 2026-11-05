# Binding-Site Detection

**Summary**: Binding-site detection identifies the region of a protein surface where a ligand is likely to bind, either from experimental co-crystal information or from computational pocket-finding and druggability analysis.
**Sources**: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Molecular dynamics simulations and drug discovery.md; Basic docking — Autodock Vina 1.2.0 documentation.md; autodock-vina-readthedocs-io-en-latest.pdf; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Molecular dynamics simulation for all.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Best Practices for Foundations in Molecular Simulations [Article.pdf
**Last updated**: 2026-11-02

---

## Why the site matters

Docking is only meaningful inside a defined search region (source: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). The binding site definition dictates the grid dimensions, the residues that may be treated with explicit conformational flexibility, and the reference coordinate frame used to judge whether a sampled pose is correct (source: Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md).

When a co-crystal structure exists, the ligand position gives the site directly. In the AutoDock Vina tutorial, the imatinib site in c-Abl is defined by the crystallographic ligand STI, and the box is centered at (15.190, 53.903, 16.917) Å with 20 Å sides (source: Basic docking — Autodock Vina 1.2.0 documentation.md). When multiple subpockets exist, defining excessive grid volumes dilutes sampling entropy, whereas overly restrictive boxes artificially truncate extended ligand conformations (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). See [[benchmark-drug-targets]].

## Sources of binding-site information

### Experimental structures

X-ray crystallography, NMR, cryo-EM, and Laue diffraction define orthosteric sites and can reveal allosteric sites occupied by known modulators. Resolution should preferably be better than 2.5 Å for quantitative work (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### Homology models and predicted structures

When no experimental structure exists, homology modelling, threading, or ab initio prediction can supply a model; AlphaFold2 and RoseTTAFold, along with all-atom co-folding foundation models like AlphaFold 3 and Chai-1, have greatly expanded coverage (source: Protein–Ligand Docking in the Machine-Learning Era.md; Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). Site detection on apo models carries additional uncertainty from side-chain and loop placement, which co-folding models aim to resolve by predicting the ligand-bound pocket directly (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### Blind docking

Blind docking lets the search explore the entire protein surface and returns favorable cavities de novo. It is computationally expensive because the search space grows with protein volume, and the top-scoring site may be an artifact of the scoring function rather than a real pocket (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Computational pocket predictors

Pocket predictors characterize the spatial composition of amino acids or use chemical probes to find favorable hot spots (source: Protein–Ligand Docking in the Machine-Learning Era.md). Commonly used tools include:

- **AlphaSpace** — maps concave surface pockets and their chemical character;
- **FTMap** — floods the surface with small organic probes and identifies regions where probes congregate;
- **MDpocket** — tracks pocket volume and openness across MD trajectories;
- **Fpocket** — geometry-based detection of alpha spheres;
- **SiteMap** — identifies and scores druggable sites by size, enclosure, hydrophobicity, and hydrogen-bonding potential.

These tools are also the entry point to allosteric and cryptic pockets, which may not be visible in any single crystal structure (source: Protein–Ligand Docking in the Machine-Learning Era.md).

### FTMap and Mixed-Solvent MD in practice

Ivetac and McCammon simulated the human β1 and β2 adrenergic receptors, extracted multiple conformations, and flooded the protein surface with small organic probes using FTMAP. Regions where probes consistently congregated across multiple structures were nominated as allosteric sites. Five potential sites were found, some absent from all existing crystal structures (source: Molecular dynamics simulations and drug discovery.md).

Furthermore, [[mixed-solvent-molecular-dynamics]] (MixMD, SILCS, SWISH, CrypticScout) advances probe mapping into full explicit-solvent dynamic simulations (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). Small organic probes (e.g., 5% v/v isopropanol, benzene, phenol, acetonitrile) compete with water across multiple independent replicates, inducing the opening of cryptic pockets and mapping 3D Grid Free Energies (GFE) $\Delta G_{\text{probe}}(\mathbf{r}) = -k_B T \ln [\rho(\mathbf{r})/\rho_{\text{bulk}}]$ to locate high-affinity hotspots for solvent-site-based docking (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). See [[mixed-solvent-molecular-dynamics]].

## Druggability

A cavity is not automatically a drug site. Useful criteria include:

- **enclosure and volume** — a pocket buried enough to exclude bulk solvent and large enough for a drug-like molecule;
- **hydrophobicity and aromatic content** — favorable for buried apolar surface;
- **hydrogen-bonding capacity** — polar groups at the pocket rim and base;
- **conservation** — conserved cavities are more likely functional;
- **competition with solvent** — waters that are hard to displace reduce attractiveness.

The **screening performance index** operationalizes part of this for VS: five calculated terms describe how well a set of structures performs on known actives. Co-crystal structures with large ligands generally score well, but the method needs labelled data that may not exist for a novel target (source: Protein–Ligand Docking in the Machine-Learning Era.md).

## Orthosteric versus allosteric and cryptic sites

- **Orthosteric site** — the endogenous ligand site. It is usually the best characterized but is prone to resistance mutations, which motivates the search for alternatives (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Allosteric site** — a distinct site whose occupation modulates activity, often by changing the population of active and inactive conformations.
- **Cryptic site** — a pocket absent or occluded in the available structures that opens transiently in solution. Cryptic sites are a major reason why static structures underestimate druggable surface (source: Molecular dynamics simulations and drug discovery.md).

See [[cryptic-and-allosteric-pockets]].

## Validation of a predicted site

A predicted site should be tested before a full screen (source: Protein–Ligand Docking in the Machine-Learning Era.md; Molecular dynamics simulations and drug discovery.md):

1. Check whether the pocket is open in multiple structures or MD frames (source: Molecular dynamics simulations and drug discovery.md).
2. Dock known actives and decoys; a site that cannot enrich known actives is suspect (source: Protein–Ligand Docking in the Machine-Learning Era.md).
3. Test whether the pocket can accommodate the intended ligand size without severe strain (source: Protein–Ligand Docking in the Machine-Learning Era.md).
4. Mutate key residues experimentally or computationally and check the effect (source: Protein–Ligand Docking in the Machine-Learning Era.md).
5. For a cryptic site, confirm with site-directed mutagenesis, covalent probes, or crystallography under different conditions (source: Molecular dynamics simulations and drug discovery.md).

## Common failure modes

- **Box too small** — part of the ligand cannot fit; Vina recommends the smallest box that still contains the whole site and warns above 27,000 Å³ (source: autodock-vina-readthedocs-io-en-latest.pdf).
- **Box too large** — search becomes inefficient and exhaustiveness must be increased.
- **Wrong protomer or missing cofactor** — the site is distorted.
- **Crystal packing artifact** — a pocket is opened or closed by lattice contacts, which MD can relax (source: Molecular dynamics simulation for all.md).
- **Scoring-driven false pocket** — a hydrophobic surface patch scores well without being druggable.

## Related pages
- [[virtual-screening]]
- [[protein-ligand-docking]]
- [[cryptic-and-allosteric-pockets]]
- [[mixed-solvent-molecular-dynamics]]
- [[gromacs-protein-ligand-protocol]]
- [[simulation-and-benchmark-best-practices]]
- [[receptor-preparation]]
- [[molecular-dynamics-simulation]]
- [[water-and-solvation]]
- [[scoring-functions]]
- [[benchmark-drug-targets]]
