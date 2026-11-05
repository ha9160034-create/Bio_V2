# Cryptic and Allosteric Pockets

**Summary**: Cryptic pockets are binding cavities that are absent or closed in available experimental structures but open transiently in solution; allosteric pockets are regulatory sites distinct from the orthosteric site. Both expand the druggable surface and both are best discovered with dynamics rather than static structures.
**Sources**: Molecular dynamics simulations and drug discovery.md; Molecular dynamics simulation for all.md; Protein–Ligand Docking in the Machine-Learning Era.md; Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf
**Last updated**: 2026-11-02

---

## Definitions

A **cryptic binding site** is a pocket that is not apparent in the available crystal or NMR structures — it may be fully closed, occluded by side chains, or simply absent from the conformational ensemble that crystallized (source: Molecular dynamics simulations and drug discovery.md).

An **allosteric site** is a regulatory site distinct from the orthosteric (endogenous ligand) site. Binding there changes the conformational equilibrium or dynamics of the protein and thereby modulates activity, often without competing directly with the substrate (source: Molecular dynamics simulations and drug discovery.md).

The categories overlap: many allosteric sites are also cryptic, because the active receptor conformation is transiently populated and not captured by a single structure.

## Why static structures hide them

X-ray crystallography and cryo-EM return an average structure, often stabilized by lattice contacts or by a bound ligand that selects one conformation. A pocket that opens only in a low-population state, or only after a loop moves, may be invisible. Molecular dynamics is well suited to finding such sites because it samples the conformational ensemble in solution (source: Molecular dynamics simulations and drug discovery.md; Molecular dynamics simulation for all.md).

## Case study: HIV integrase and raltegravir

This is the canonical cryptic-pocket success. In 2004, Schames and coworkers ran MD on HIV integrase, a target that had seemed intractable for structure-based design. The simulation revealed a previously unidentified trench that was not evident in any crystal structure. X-ray crystallography later showed that known inhibitors bind in that trench, as predicted. The result led to new experimental work at Merck, and further development produced **raltegravir**, the first FDA-approved HIV integrase inhibitor (source: Molecular dynamics simulations and drug discovery.md).

The lesson is that the simulation did not merely refine a known site; it generated a structural hypothesis that changed a drug program.

## Case study: β-adrenergic receptor allosteric sites

Ivetac and McCammon simulated the human β1 and β2 adrenergic receptors. They extracted multiple conformations and flooded the protein surface with small organic probes using FTMAP. Regions where probes consistently congregated across multiple structures were nominated as allosteric sites. Five potential sites were identified, some absent from all existing crystal structures (source: Molecular dynamics simulations and drug discovery.md).

## Computational routes to cryptic sites

### 1. Long unbiased MD

Multi-microsecond simulations can capture a pocket opening spontaneously. The D. E. Shaw group's Anton simulations of Src kinase showed PP1 docking at alternative pockets in addition to the orthosteric site, demonstrating that long unbiased runs can reveal allosteric sites not apparent crystallographically (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### 2. Pocket tracking on trajectories

MDpocket, FTMap, and related tools can be run on trajectory frames to detect cavities that appear and disappear. The key output is not a single structure but a **pocket occupancy** — the fraction of frames in which a cavity of adequate volume is open (source: Protein–Ligand Docking in the Machine-Learning Era.md; Molecular dynamics simulations and drug discovery.md).

### 3. Enhanced sampling

Accelerated MD, metadynamics, scaled MD, and weighted ensemble methods increase the probability of visiting high-free-energy states where pockets open (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[enhanced-sampling]].

### 4. Markov state models

MSMs built from many short trajectories identify metastable conformational states and their populations. Different apo states can expose different pockets; for β-trypsin, six apo conformations visited on tens-of-microseconds timescales exposed distinct surface pockets, and these conformations were recognizable in crystal structures of wild-type and mutant protein (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

### 5. Mixed solvent and probe methods (MixMD / SILCS / SWISH)

Small organic probes in explicit cosolvent MD simulations actively induce cryptic pocket opening while mapping surface hotspots (source: Molecular dynamics simulations and drug discovery.md; Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). In [[mixed-solvent-molecular-dynamics]] (MixMD, SILCS, SWISH, CrypticScout), organic probes (typically 5% v/v isopropanol, benzene, phenol, or acetonitrile) compete with water across multiple independent replicates; probe insertion lowers the free energy barrier for loop displacement, exposing cryptic pockets (such as the sterol-binding pocket of Niemann-Pick Type C2 protein) and mapping 3D Grid Free Energies (GFE) $\Delta G_{\text{probe}}(\mathbf{r}) = -k_B T \ln [\rho(\mathbf{r})/\rho_{\text{bulk}}]$ to guide solvent-site-based docking (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). See [[mixed-solvent-molecular-dynamics]].

## Using a cryptic pocket in docking

Because rigid-receptor docking protocols cannot accommodate sidechain rearrangements or loop movements on the fly, static docking against closed apo crystal structures frequently produces catastrophic false negatives (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md). Once a cryptic or allosteric conformation is sampled, it can be exploited in structure-based pipelines:

- **ensemble docking** — add the cryptic conformation to the receptor ensemble and dock the library against all conformers, mitigating rigid-receptor bias across conformational states (source: Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md; Molecular dynamics simulations and drug discovery.md);
- **the relaxed complex scheme (RCS)** — dock each ligand into many MD-derived conformations and rank by a spectrum statistic such as the average score. RCS has identified inhibitors of FKBP, HIV integrase, *Trypanosoma brucei* RNA editing ligase 1, *T. brucei* GalE, *T. brucei* FPPS, and *Mycobacterium tuberculosis* dTDP-6-deoxy-L-lyxo-4-hexulose synthase; in two projects the inhibitors were active against the whole parasite (source: Molecular dynamics simulations and drug discovery.md);
- **free-energy calculation** — compute absolute binding free energy to the cryptic conformation to test whether binding is actually favorable (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[benchmark-drug-targets]].

## Practical cautions

- **A pocket opening is not proof of druggability.** The open state may be too rare, too polar, or too shallow.
- **Force-field artifacts.** Pocket opening can be exaggerated by imperfect force fields or by the absence of a membrane or binding partner (source: Molecular dynamics simulations and drug discovery.md).
- **Sampling is incomplete.** A simulation that fails to open a pocket does not prove the pocket does not exist; results are incomplete rather than wrong (source: Molecular dynamics simulations and drug discovery.md).
- **Ensemble cost and false positives.** Adding many receptor conformers increases both computational cost and false-positive rate in VS (source: Protein–Ligand Docking in the Machine-Learning Era.md).
- **Validation is essential.** Site-directed mutagenesis, covalent probes, NMR relaxation, or new crystal forms should test the predicted site.

## Allostery as a thermodynamic mechanism

MD simulations of the M2 muscarinic acetylcholine receptor showed how allosteric modulators change the conformation of the orthosteric site: the orthosteric ligand N-methyl scopolamine favors an enlarged allosteric site; the positive modulator alcuronium requires the larger site, whereas the negative modulator does not. A designed modulator predicted to need a larger pocket had fourfold less negative cooperativity, even though its affinity was higher — a direct experimental validation of a simulation-derived allosteric mechanism (source: Molecular dynamics simulation for all.md).

## Related pages
- [[binding-site-detection]]
- [[mixed-solvent-molecular-dynamics]]
- [[gromacs-protein-ligand-protocol]]
- [[simulation-and-benchmark-best-practices]]
- [[molecular-dynamics-simulation]]
- [[enhanced-sampling]]
- [[virtual-screening]]
- [[structure-based-drug-design]]
- [[binding-free-energy]]
- [[benchmark-drug-targets]]
- [[durrant-2011-md-drug-discovery]]
- [[markov-state-models]]
