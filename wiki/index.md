# Wiki Index

**Summary**: Table of contents for a knowledge base on protein–ligand docking, molecular dynamics, binding thermodynamics and kinetics, and computer-aided drug discovery.
**Last updated**: 2026-11-02

---

## Start here

- [[protein-ligand-docking]] — Core definition of docking, the pose-versus-affinity problem, and the sampling/scoring pipeline.
- [[structure-based-drug-design]] — How docking fits into hit finding, lead optimization, and preclinical development.
- [[molecular-dynamics-simulation]] — What MD computes, how it is run, what it can answer, and where it fails.

## Docking

- [[molecular-docking-algorithms]] — Search strategies: incremental construction, conformer libraries, Monte Carlo, genetic algorithms, and flexible-receptor methods.
- [[autodock-vina]] — The Vina docking engine: scoring function, search algorithm, practical settings, and the 1.2.0 feature set.
- [[binding-site-detection]] — Finding and validating orthosteric, allosteric, and cryptic pockets; box definition and druggability.
- [[cryptic-and-allosteric-pockets]] — Pockets absent from static structures, the raltegravir story, and MD-based discovery.
- [[ligand-preparation-and-protonation]] — SMILES/SDF to docking-ready 3D structures: tautomers, protonation, charges, and PDBQT.
- [[receptor-preparation]] — Choosing, cleaning, protonating, and formatting a receptor; waters, metals, cofactors, and flexible side chains.
- [[docking-failure-modes]] — Physical, algorithmic, and chemical causes of failure in docking and virtual screening: pocket rigidity, MW bias, ligand strain, metal coordination, and false-positive artifacts.
- [[benchmark-drug-targets]] — Reference therapeutic case studies across CADD: HIV-1 protease/integrase, muscarinic M3 receptor, kinase families (c-Abl, Src, p38, CDK8, MAPK2), and COMT metalloenzyme coordination.

## Scoring and machine learning

- [[scoring-functions]] — Force-field, empirical, and knowledge-based scoring; what each can and cannot predict.
- [[machine-learning-scoring-functions]] — From RF-Score to graph neural networks; delta-ML, interaction features, and generalization limits.
- [[deep-learning-docking]] — Generative diffusion, regression, and co-folding docking models (DiffDock, AF3, Chai-1, Boltz-1), PoseBusters validity, dataset memorization, and cross-docking collapse.
- [[benchmarking-datasets]] — PDBbind, CASF, DUD-E, PoseBusters, DockGen, Astex Diverse, MUV, LIT-PCBA, CSAR, D3R, and CACHE.

## Virtual screening and design

- [[virtual-screening]] — The five-step structure-based VS workflow, ultra-large libraries, enrichment metrics, chemical artifacts (PAINS, aggregators), and multi-tier rescoring cascades.
- [[chemical-screening-artifacts]] — Promiscuous false-positive mechanisms (PAINS electrophiles, colloidal aggregation), detergent counter-screens (Triton X-100), and in silico/in vitro artifact remediation.

## Molecular dynamics and free energy

- [[force-fields]] — Bonded and non-bonded terms, AMBER/CHARMM/OPLS families, water models, polarization, and QM/MM.
- [[gromacs-protein-ligand-protocol]] — Complete step-by-step GROMACS simulation protocol: topology generation, CGenFF/ACPYPE ligand parameterization, position restraints, merged thermostat groups, CHARMM36 force-switch settings, and interaction energy decomposition.
- [[open-force-field-and-ligand-parameterization]] — Modern small-molecule parameterization: direct chemical perception (SMIRNOFF), OpenFF Sage 2.2.0, continuous neural network force fields (Garnet), and comparison with legacy atom-typing workflows (GAFF, CGenFF, OPLS).
- [[mixed-solvent-molecular-dynamics]] — Explicit cosolvent MD (MixMD, SILCS): organic probe sweeps (5% v/v), 3D Grid Free Energy (GFE) hotspot mapping, cryptic pocket induction, and integration into solvent-site-based docking.
- [[simulation-and-benchmark-best-practices]] — Consolidated LiveCoMS best-practice checklists for foundational MD simulation validation (integrators, thermostats, barostats, PBC) and Hahn et al. standards for constructing, preparing, and evaluating protein–ligand binding affinity benchmarks.
- [[hydrogen-mass-repartitioning]] — Biophysical protocol and numerical basis of HMR, transferring mass to hydrogens to enable stable 4.0 fs timesteps in MD without altering thermodynamic ensembles.
- [[enhanced-sampling]] — CV-based vs non-CV methods, Steered MD, ABMD, parallel tempering, scaled MD, τ-RAMD, metadynamics, GaMD, SuMD, and weighted ensemble.
- [[metadynamics]] — History-dependent repulsive bias deposition along collective variables, Well-Tempered convergence, and Infrequent Metadynamics (InMetaD) for drug unbinding kinetics.
- [[umbrella-sampling]] — Stratified equilibrium sampling along reaction coordinates using harmonic bias windows, WHAM/MBAR unbiasing, and standard-state volume corrections ($V_0 = 1661\text{ \AA}^3$).
- [[markov-state-models]] — Master-equation kinetic modeling from distributed MD ensembles: TICA dimensionality reduction, microstate clustering, Chapman–Kolmogorov validation, and transition path flux.
- [[binding-free-energy]] — The computational binding hierarchy from docking to endpoint continuum solvent and alchemical perturbation.
- [[mm-pbsa-and-mm-gbsa]] — Endpoint binding free energy methods: single- vs three-trajectory MD, PB vs GB continuum electrostatics, dielectric constant tuning, and virtual screening rescoring.
- [[alchemical-free-energy-methods]] — Thermodynamic cycles, FEP, TI, BAR, MBAR, soft-core potentials, charge decoupling, buried waters, active learning, and neural network potentials.
- [[water-and-solvation]] — The hydrophobic effect, bridging waters, desolvation barriers, and how docking and MD treat solvent.

## Thermodynamics and kinetics

- [[binding-thermodynamics]] — ΔG, ΔH, ΔS, enthalpy–entropy compensation, ITC, SPR, fluorescence polarization, and heat capacity.
- [[experimental-binding-assays]] — Principles, operational instrumentation, and deconvolution of thermodynamic and kinetic parameters via ITC and SPR, with methods to detect experimental artifacts.
- [[binding-kinetics]] — k_on, k_off, residence time, clinical paradigms (tiotropium vs ipratropium), and computational estimates from unbiased MD, InMetaD, LiGaMD, SuMD, and τ-RAMD.
- [[protein-ligand-interactions]] — Hydrogen bonds, ionic pairs, van der Waals contacts, hydrophobic burial, halogen bonds, and metal coordination.
- [[molecular-recognition-models]] — Lock-and-key, induced fit, conformational selection, and extended conformational selection.

## Source summaries

- [[ma-2024-docking-strategies]] — Zou-group review of docking physics, sampling, scoring, programs, inverse VS, and the CELPP 80.4%/46.3% gap.
- [[huang-2010-scoring-functions]] — Foundational review of the three scoring-function families, consensus scoring, and evaluation criteria.
- [[yang-2022-ml-docking]] — Review of ML scoring functions, datasets, the VS workflow, and the ΔLin_F9 XGB LIT-PCBA benchmark.
- [[durrant-2011-md-drug-discovery]] — MD for cryptic sites, the relaxed complex scheme, and alchemical free energy; the raltegravir case.
- [[hollingsworth-2018-md-for-all]] — Accessible review of MD practice: time steps, force fields, software, applications, and limitations.
- [[du-2016-protein-ligand-interactions]] — Thermodynamics, kinetics, recognition models, ITC/SPR/FP, and the HIV protease inhibitor examples.
- [[decherchi-2020-binding-thermodynamics-kinetics]] — Comprehensive review of enhanced sampling, free-energy methods, kinetics, and practical guidelines.
- [[vina-docs-basic-docking]] — The Vina 1.2.0 manual and tutorials: 1IEP basic docking, flexible receptor, multiple ligands, zinc, water, macrocycles, and Python.
- [[biochemistry-free-for-all]] — OER biochemistry textbook; front matter and scope only, with detailed chapters not yet ingested.
- [[gimeno-2019-virtual-screening-pitfalls]] — Critical review of virtual screening methodologies, benchmark decoy biases, PAINS, colloidal aggregation, and orthogonal assay validation.
- [[salmaso-2018-bridging-docking-md]] — Review bridging molecular docking to MD: molecular recognition models, ensemble docking, enhanced sampling, and Supervised MD (SuMD).
- [[xu-2022-docking-failures-analysis]] — DOCK 3.7 vs AutoDock Vina benchmark on DUD-E: 23x speed difference, Vina's empirical MW bias, 311 irregular torsions, and COMT metalloenzyme collapse.
- [[genheden-2015-mm-pbsa-gbsa]] — Authoritative review of MM/PBSA and MM/GBSA: 1-trajectory vs 3-trajectory, internal dielectric calibration, non-polar terms, and the benefit of omitting normal mode entropy.
- [[poli-2020-mm-pbsa-virtual-screening]] — Application of MM-PBSA in virtual screening: energy minimization vs MD protocols, 1% Enrichment Factor improvements, and dielectric calibration.
- [[bernardi-2015-enhanced-sampling-md]] — Schulten group review of enhanced sampling in biological MD: Replica Exchange (REMD), Metadynamics, and Generalized Simulated Annealing (GSA).
- [[lazim-2020-md-enhanced-sampling-advances]] — Classification of enhanced sampling into CV-based and non-CV-based (GaMD), Hydrogen Mass Repartitioning (4 fs timestep), and GPCR/kinase applications.
- [[york-2023-modern-alchemical-free-energy]] — Detailed physical explanation of alchemical free energy: TI, FEP, BAR, MBAR, soft-core potentials, staged electrostatics, and error diagnostics.
- [[muegge-2023-alchemical-free-energy-advances]] — Industry perspective on alchemical free energy: RBFE vs ABFE, fragment elaboration, 30–50% identity homology models, buried waters with GCMC, and covalent inhibitors.
- [[van-pinxteren-2025-ml-fep-workflows]] — Review of machine learning integration into FEP: active learning with up to 20-fold screening speedups, AIM-Net2 ML potentials, and DL holo structures.
- [[morehead-2026-evaluating-dl-docking]] — Nature Machine Intelligence evaluation of deep learning docking (DiffDock, Boltz-1, Chai-1, AF3) vs P2Rank-Vina on PoseBusters PB-valid compliance and clash rates.
- [[li-2025-decoding-dl-docking-limits]] — Chemical Science benchmark uncovering 0% physical validity in KarmaDock, ligand memorization on DockGen, cross-docking collapse, and virtual screening failure.
- [[ahmad-2022-enhanced-sampling-binding-kinetics]] — Review of enhanced-sampling MD for unbinding kinetics ($k_{\text{off}}$) and residence time: InMetaD, OPES-flooding, $\tau$-RAMD, LiGaMD, and Markov State Models.
- [[adediwura-2024-free-energy-kinetics-drug-discovery]] — Review from the Miao lab unifying alchemical free energy, MM-PBSA/GBSA, LiGaMD kinetics, and clinical residence time paradigms (tiotropium vs ipratropium).
- [[lemkul-gromacs-protein-ligand-tutorial]] — Justin Lemkul's authoritative tutorial detailing the complete simulation workflow for a protein–ligand complex (T4 lysozyme + JZ4) in GROMACS using CHARMM36 and CGenFF.
- [[lemkul-2018-gromacs-tutorials]] — Lemkul's peer-reviewed LiveCoMS article presenting a suite of seven standardized GROMACS tutorials spanning aqueous proteins, membranes, umbrella sampling, biphasic systems, ligands, and alchemical free energy.
- [[braun-2019-foundations-molecular-simulation]] — LiveCoMS foundational best practices guide for molecular simulations: checklist, statistical mechanical ensembles, thermostat/barostat audit, and PBC artifacts.
- [[hahn-2022-binding-affinity-benchmarks]] — LiveCoMS best practices guide for constructing, preparing, and evaluating protein–ligand binding affinity benchmarks (PLBenchmarks and arsenic).
- [[wang-2024-open-force-field-initiative]] — Status report on the Open Force Field Initiative: SMIRNOFF direct chemical perception, ForceBalance optimization, and the evolution of open-source biomolecular force fields from Parsley to Sage 2.2.0.
- [[blanco-gonzalez-2024-garnet-force-field]] — MRC LMB research presenting Garnet: an end-to-end graph neural network parameterizing classical force fields from scratch for proteins and small molecules using double exponential nonbonded potentials.
- [[poziemski-2026-md-time-series-descriptors]] — Large-scale benchmark across 2,752 complexes evaluating machine learning models on MD time-series descriptors, demonstrating superior target and scaffold out-of-distribution generalizability.

## Reading paths

**Beginner: Foundations of Docking and Screening**: [[protein-ligand-docking]] → [[molecular-docking-algorithms]] → [[scoring-functions]] → [[autodock-vina]] → [[vina-docs-basic-docking]] → [[gimeno-2019-virtual-screening-pitfalls]]

**Intermediate: Practical Workflows, Failures & Rescoring**: [[virtual-screening]] → [[chemical-screening-artifacts]] → [[docking-failure-modes]] → [[benchmark-drug-targets]] → [[xu-2022-docking-failures-analysis]] → [[salmaso-2018-bridging-docking-md]] → [[mm-pbsa-and-mm-gbsa]] → [[genheden-2015-mm-pbsa-gbsa]] → [[poli-2020-mm-pbsa-virtual-screening]]

**Advanced: Rigorous Free Energy, Kinetics & Deep Learning**: [[alchemical-free-energy-methods]] → [[york-2023-modern-alchemical-free-energy]] → [[muegge-2023-alchemical-free-energy-advances]] → [[van-pinxteren-2025-ml-fep-workflows]] → [[simulation-and-benchmark-best-practices]] → [[hahn-2022-binding-affinity-benchmarks]] → [[deep-learning-docking]] → [[morehead-2026-evaluating-dl-docking]] → [[li-2025-decoding-dl-docking-limits]] → [[binding-kinetics]] → [[experimental-binding-assays]] → [[markov-state-models]] → [[metadynamics]] → [[ahmad-2022-enhanced-sampling-binding-kinetics]] → [[adediwura-2024-free-energy-kinetics-drug-discovery]]

**MD and Enhanced Sampling Specialist**: [[molecular-dynamics-simulation]] → [[force-fields]] → [[open-force-field-and-ligand-parameterization]] → [[wang-2024-open-force-field-initiative]] → [[blanco-gonzalez-2024-garnet-force-field]] → [[gromacs-protein-ligand-protocol]] → [[lemkul-gromacs-protein-ligand-tutorial]] → [[lemkul-2018-gromacs-tutorials]] → [[simulation-and-benchmark-best-practices]] → [[braun-2019-foundations-molecular-simulation]] → [[mixed-solvent-molecular-dynamics]] → [[hydrogen-mass-repartitioning]] → [[enhanced-sampling]] → [[umbrella-sampling]] → [[metadynamics]] → [[bernardi-2015-enhanced-sampling-md]] → [[lazim-2020-md-enhanced-sampling-advances]] → [[binding-free-energy]]
