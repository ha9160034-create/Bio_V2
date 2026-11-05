# Wiki Log

Append-only record of wiki operations. Newest entries at the bottom.

---

## 2026-11-01 — Initial full ingestion of the knowledge base

**Sources read**

All ten files in `row/` (the repository has no `raw/` folder; `row/` was treated as the immutable source collection, and its filenames are used verbatim in citations):

1. `Docking strategies for predicting protein-ligand interactions and their application to structure-based drug design.md`
2. `Scoring functions and their evaluation methods for protein-ligand docking_ recent advances and future directions.md`
3. `Protein–Ligand Docking in the Machine-Learning Era.md` (en dash in `Protein–Ligand`)
4. `Molecular dynamics simulations and drug discovery.md`
5. `Molecular dynamics simulation for all.md`
6. `Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md` (en dash in `Protein–Ligand`)
7. `Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md`
8. `Basic docking — Autodock Vina 1.2.0 documentation.md` (em dash)
9. `autodock-vina-readthedocs-io-en-latest.pdf`
10. `Biochemistry Free For All 1.3_compressed.pdf`

**Processing**

- Markdown sources were read in full. Non-ASCII characters in filenames were preserved: U+2013 in the Yang and Du filenames, U+2014 in the Vina basic-docking filename.
- The Vina PDF was extracted in full with Python `pypdf` 6.21.4 (44 pages, ~89,000 characters).
- The biochemistry PDF is 3,648 pages. Only the front matter, table of contents, and roughly the first 40 pages were extracted (~58,000 characters). The detailed chapters were not ingested.

**Pages created — concepts (21)**

- `protein-ligand-docking.md`
- `molecular-docking-algorithms.md`
- `scoring-functions.md`
- `machine-learning-scoring-functions.md`
- `autodock-vina.md`
- `structure-based-drug-design.md`
- `virtual-screening.md`
- `binding-site-detection.md`
- `cryptic-and-allosteric-pockets.md`
- `protein-ligand-interactions.md`
- `molecular-recognition-models.md`
- `molecular-dynamics-simulation.md`
- `force-fields.md`
- `enhanced-sampling.md`
- `binding-free-energy.md`
- `binding-thermodynamics.md`
- `binding-kinetics.md`
- `water-and-solvation.md`
- `ligand-preparation-and-protonation.md`
- `receptor-preparation.md`
- `benchmarking-datasets.md`

**Pages created — source summaries (9)**

- `ma-2024-docking-strategies.md`
- `huang-2010-scoring-functions.md`
- `yang-2022-ml-docking.md`
- `durrant-2011-md-drug-discovery.md`
- `hollingsworth-2018-md-for-all.md`
- `du-2016-protein-ligand-interactions.md`
- `decherchi-2020-binding-thermodynamics-kinetics.md`
- `vina-docs-basic-docking.md`
- `biochemistry-free-for-all.md`

**Pages created — navigation**

- `index.md` — grouped table of contents with one-line descriptions and suggested reading paths.
- `log.md` — this file.

**Decisions and caveats**

- `row/` was used instead of the `raw/` directory named in `AGENTS.md` because no `raw/` exists. `row/` was treated as read-only and was not modified.
- The 3,648-page biochemistry textbook was only partially read. Its page explicitly flags this limitation and marks detailed content as needing verification.
- Several source files contain typographic dashes in their names. Citations in wiki pages use the exact filename characters where the page was written from scratch; earlier-created pages may cite a plain hyphen in `Protein–Ligand Docking in the Machine-Learning Era.md`, `Insights into Protein–Ligand Interactions_...md`, or `Basic docking - ...md`. A normalization pass was planned to make all citations match the stored filenames exactly.
- Two links in `protein-ligand-docking.md` — `ensemble-docking` and `soft-docking` — had no dedicated pages. They are covered inside `virtual-screening.md` and `molecular-docking-algorithms.md`; the links were redirected rather than creating thin stub pages.
- `wiki-links` in this file is written as code, not a link, because it is a syntax example rather than a page.

**Verification**

- Every page was checked for the required header block (`# Title`, `**Summary**`, `**Sources**`, `**Last updated**`), a `---` divider, and a `## Related pages` section.
- All internal wiki links were checked against the actual filenames in `wiki/`; broken links were repaired or redirected to the pages that cover the topic.

**Next steps**

- Ingest selected biochemistry chapters (Water and Buffers, Amino Acids, Protein Structure, Catalysis, Energy) if textbook-level definitions are needed.
- Add a dedicated page for ensemble docking and soft docking if those topics grow beyond their current treatment in `virtual-screening.md`.
- Add prospective benchmark results as they become available (CACHE, D3R).

---

## 2026-11-01 — Citation normalization and verification pass

**Changes**

- Normalized source filenames in citations to match the stored files exactly: U+2013 en dash in `Protein–Ligand Docking in the Machine-Learning Era.md` and `Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md`; U+2014 em dash in `Basic docking — Autodock Vina 1.2.0 documentation.md`.
- Added `autodock-vina-readthedocs-io-en-latest.pdf` as the cited source where facts came from the full manual (scoring, search, box warnings, flexible docking, multiple ligands, zinc, hydrated docking, macrocycles, Python API) rather than the basic-docking page.
- Added every cited source to the corresponding page's `**Sources**` line; the consistency check now reports zero pages citing an unlisted source.
- Fixed the `**Last updated**` header in `water-and-solvation.md`.
- Removed unresolved `ensemble-docking` and `soft-docking` links from `protein-ligand-docking.md`; the topics are covered in `molecular-docking-algorithms.md` and `virtual-screening.md`.

**Final verification**

- 32 pages total: 21 concept pages, 9 source summaries, `index.md`, and `log.md`.
- 385 internal wiki links; **0 broken**.
- All 30 content pages pass the header, divider, and `## Related pages` checks.
- No content page is orphaned (every page has at least one inbound link).
- The biochemistry textbook remains partially ingested, as flagged on its source page.

---

## 2026-11-01 — Final accuracy polish

**Changes**

- Corrected the relative FEP/TI accuracy entry in `binding-free-energy.md` from "~1 kcal/mol; 81.2% within 1.5 kcal/mol" to "≤1.5 kcal/mol in 81.2% of cases", matching Decherchi & Cavalli exactly.
- Reworded a verification bullet in this log so placeholder link names appear as plain code rather than as resolvable wiki links.

**Final verification**

- 32 pages; 385 internal links; 0 broken; no orphaned content pages.
- All 30 content pages pass the header, divider, and `## Related pages` checks.
- `index.md` links every content page.


---

## 2026-11-01 — Systematic ingestion of 15 papers across beginner, intermediate, and advanced levels

**Sources read and ingested**

All 15 new PDF sources in `row/` were extracted and processed systematically:

1. `The Light and Dark Sides of Virtual Screening.pdf` (Gimeno et al., 2019, IJMS) — Beginner/Foundational.
2. `Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf` (Salmaso & Moro, 2018, Front. Pharmacol.) — Beginner/Foundational.
3. `Systematic Investigation of Docking Failures in Large-Scale Structure-Based Virtual Screening.pdf` (Xu et al., 2022, ACS Omega) — Intermediate.
4. `The MMPBSA and MMGBSA methods to estimate ligand-binding affinities.pdf` (Genheden & Ryde, 2015, Expert Opin. Drug Discov.) — Intermediate.
5. `Application of MM-PBSA Methods in Virtual Screening.pdf` (Poli et al., 2020, Molecules) — Intermediate.
6. `Enhanced Sampling Techniques in Molecular Dynamics Simulations of Biological Systems.pdf` (Bernardi, Melo & Schulten, 2015, BBA) — Intermediate/Advanced.
7. `Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf` (Lazim, Suh & Choi, 2020, IJMS) — Intermediate/Advanced.
8. `Modern Alchemical Free Energy Methods for Drug Discovery Explained.pdf` (York, 2023, ACS Phys. Chem Au) — Advanced.
9. `Recent Advances in Alchemical Binding Free Energy Calculations for Drug Discovery.pdf` (Muegge & Hu, 2023, ACS Med. Chem. Lett.) — Advanced.
10. `Integrating Machine Learning into Free Energy Perturbation Workflows.pdf` (van Pinxteren & Jespers, 2025, JCIM) — Advanced.
11. `Assessing the potential of deep learning for protein–ligand docking.pdf` (Morehead et al., 2026, Nature Machine Intelligence) — Advanced.
12. `Decoding the limits of deep learning in molecular docking for drug discovery.pdf` (Li et al., 2025, Chemical Science) — Advanced.
13. `Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf` (Ahmad et al., 2022, Front. Mol. Biosci.) — Advanced.
14. `Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf` (Adediwura et al., 2024, Expert Opin. Drug Discov.) — Advanced.
15. `Scoring functions and their evaluation methods for protein-ligand docking recent advances and future directions.pdf` (Huang, Grinter & Zou, 2010, PCCP) — Added to existing summary `huang-2010-scoring-functions.md`.

**Pages created — source summaries (14)**

- `gimeno-2019-virtual-screening-pitfalls.md`
- `salmaso-2018-bridging-docking-md.md`
- `xu-2022-docking-failures-analysis.md`
- `genheden-2015-mm-pbsa-gbsa.md`
- `poli-2020-mm-pbsa-virtual-screening.md`
- `bernardi-2015-enhanced-sampling-md.md`
- `lazim-2020-md-enhanced-sampling-advances.md`
- `york-2023-modern-alchemical-free-energy.md`
- `muegge-2023-alchemical-free-energy-advances.md`
- `van-pinxteren-2025-ml-fep-workflows.md`
- `morehead-2026-evaluating-dl-docking.md`
- `li-2025-decoding-dl-docking-limits.md`
- `ahmad-2022-enhanced-sampling-binding-kinetics.md`
- `adediwura-2024-free-energy-kinetics-drug-discovery.md`

**Pages created — concepts (4)**

- `mm-pbsa-and-mm-gbsa.md` — Intermediate endpoint binding free energy methods: single- vs three-trajectory approach, continuum electrostatics (PB vs GB), non-polar terms (SASA vs cavity/dispersion), dielectric constant calibration ($\epsilon_{\text{int}} = 1, 2, 4$), entropy calculation dilemmas (NMA vs IE vs omission), and multi-tier virtual screening rescoring protocols.
- `alchemical-free-energy-methods.md` — Advanced alchemical free energy calculations: thermodynamic cycles, statistical mechanics estimators (TI, FEP, BAR, MBAR, NEW), soft-core potential formulations, staged charge decoupling, single vs dual topology, buried water sampling via GCMC, net charge finite-size corrections, and active-learning/ML potential integration.
- `deep-learning-docking.md` — Advanced deep learning docking architectures: generative diffusion models (DiffDock, DynamicBind, DiffBindFR), regression models (EquiBind, TANKBind, FABind), all-atom co-folding foundation models (AlphaFold 3, Chai-1, Boltz-1, RFAA), PoseBusters physical validity criteria (PB-Valid), dataset memorization on PDBbind, cross-docking and apo collapse, and virtual screening failure.
- `docking-failure-modes.md` — Cross-level synthesis of docking and virtual screening failure modes: receptor rigidity and loop deformation (MAPK2 80s loop), ligand internal strain and irregular torsions (311/281 types in Vina/DOCK), empirical molecular weight scoring bias (Vina fragment penalty), metalloenzyme coordination errors (COMT $Mg^{2+}$ collapse), and chemical artifacts (PAINS, colloidal aggregators, DUD-E decoy bias).

**Concept pages updated (8)**

- `binding-free-energy.md` — Updated with modern estimators (MBAR), staged electrostatics, homology model applicability (30–50% identity matching crystal accuracy), active learning FEP (up to 20x speedup), neural network potentials (AIM-Net2), and links to new concept pages.
- `binding-kinetics.md` — Updated with clinical residence time paradigms (tiotropium vs ipratropium on muscarinic $M_3$, kinase inhibitor selectivity), Infrequent Metadynamics (InMetaD, Poisson statistics), Ligand GaMD (LiGaMD), Supervised MD (SuMD), and OPES-flooding.
- `virtual-screening.md` — Updated with DOCK 3.7 vs AutoDock Vina benchmark (23x speedup, empirical MW scoring bias, failure modes), PAINS and colloidal aggregation detection, and multi-tier MM-PBSA rescoring workflows (EF1% improvements from ~1.2 to >5.0).
- `enhanced-sampling.md` — Updated with CV-based vs non-CV-based categorization, Gaussian Accelerated MD (GaMD) and LiGaMD (2nd-order cumulant expansion reweighting), Hydrogen Mass Repartitioning (HMR, 4 fs timestep, 2-fold acceleration), Generalized Simulated Annealing (GSA), and Supervised MD (SuMD).
- `machine-learning-scoring-functions.md` — Updated with the evolution into end-to-end DL docking models, PoseBusters physical validity limits, and generalization breakdowns on novel folds (DockGen).
- `protein-ligand-docking.md` — Updated with physics-based vs deep learning comparisons, PB-valid criteria, cross-docking vs redocking gaps, and links to failure modes and alchemical methods.
- `benchmarking-datasets.md` — Updated with modern evaluation benchmarks: PoseBusters Benchmark, Astex Diverse Set, DockGen/DockGen-E, and CASP15 PLI.
- `molecular-dynamics-simulation.md` — Updated with Hydrogen Mass Repartitioning (HMR, 4.0 fs timestep) and bridging docking poses to MD relaxation and endpoint MM-PBSA rescoring.
- `autodock-vina.md` — Updated with empirical benchmark failure modes from Xu et al. (2022): MW scoring bias, lack of internal strain restrictions, and metalloenzyme charge handling.
- `huang-2010-scoring-functions.md` — Added `Scoring functions and their evaluation methods for protein-ligand docking recent advances and future directions.pdf` to sources and citations.

**Navigation and verification**

- `index.md` updated with all 18 new pages and reorganized into Beginner, Intermediate, and Advanced reading paths.
- Total wiki pages: **50 pages** (25 concept pages, 23 source summaries, `index.md`, `log.md`).
- Automated lint check passed: **0 broken links**, **0 orphaned content pages**, 100% adherence to AGENTS.md template and citation rules.

---

## 2026-11-01 — Comprehensive wiki lint remediation and concept expansion pass

**Context and scope**

In response to the full wiki audit and lint report against `AGENTS.md` guidelines, executed targeted remediations across the entire wiki: resolved header-to-citation discrepancies without deleting unreferenced sources (Finding #1), contextualized contradictions and outdated claims with modern benchmark data (Findings #4, #5, #6, #7), created 3 missing foundational concept pages (Findings #9, #10, #11), and cross-linked all pages.

**Pages created — concepts (3)**

1. `experimental-binding-assays.md` — Quantitative biophysical binding assays (ITC and SPR): operational instrumentation (adiabatic feedback vs surface plasmon sensorgrams), deconvolution of thermodynamic ($\Delta G^\circ, \Delta H^\circ, -T\Delta S^\circ, \Delta C_p$) and kinetic ($k_\text{on}, k_\text{off}, t_R$) parameters, Wiseman $c$-value optimal windows ($10 < c < 500$), mass transport limitations, and experimental artifact mitigation (detergent tests for colloidal aggregation, non-specific binding).
2. `hydrogen-mass-repartitioning.md` — Protocol and physical basis of HMR in classical MD: mass reallocation formula from heavy atoms to bound hydrogens ($\sim 3.0\text{ amu}$), lowering of harmonic vibrational frequencies ($\omega = \sqrt{k/\mu}$), holonomic constraints (LINCS/SHAKE), enabling numerically stable 4.0 fs timesteps (2-fold throughput gain) without altering configuration-space partition functions, thermodynamic ensembles, or lipid order parameters.
3. `benchmark-drug-targets.md` — Synthesized reference of recurring therapeutic case studies:
   - *HIV-1 Protease & Integrase*: Flap dynamics, structural catalytic water Wat301 displacement by cyclic ureas, and cryptic pocket discovery (relaxed complex scheme) yielding raltegravir.
   - *Muscarinic $M_3$ Receptor*: Clinical binding kinetics paradigm comparing tiotropium ($\tau \approx 34.7\text{ h}$, once-daily dosing, subtype selectivity) versus ipratropium ($\tau \approx 0.26\text{ h}$, four-times-daily dosing).
   - *Kinase Families (c-Abl, Src, p38, CDK8, MAPK2)*: DFG-in/out conformational switching (Type I vs Type II inhibition, imatinib/dasatinib/BIRB-796 kinetics), MAPK2 80s loop steric occlusion in rigid docking, and CDK8 deep pocket FEP validation.
   - *Catechol-O-methyltransferase (COMT)*: Octahedral catalytic $\text{Mg}^{2+}$ coordination failure in standard Vina ($EF_{1\%} \approx 0$) and resolution via QM/MM charge redistribution or AutoDock4Zn parameters.

**Header-to-citation discrepancy fixes (17 pages)**

Preserved all listed sources in `**Sources**:` header lines and added explicit inline `(source: ...)` citations grounded in their specific contributions:
- `alchemical-free-energy-methods.md`: Added citations for Bruce et al. (2018) and Durrant & McCammon (2011).
- `benchmarking-datasets.md`: Added citations for Ferreira et al. (2015), Gimeno et al. (2019), and Xu et al. (2022).
- `binding-free-energy.md`: Added citations for Hospital et al. (2015) and Adediwura et al. (2024).
- `binding-site-detection.md`: Added citation for Ferreira et al. (2015).
- `binding-thermodynamics.md`: Added citation for Ferreira et al. (2015).
- `cryptic-and-allosteric-pockets.md`: Added citation for Ferreira et al. (2015).
- `deep-learning-docking.md`: Added citation for Yang et al. (2022).
- `enhanced-sampling.md`: Added citations for Durrant & McCammon (2011) and Hospital et al. (2015).
- `force-fields.md`: Added citation for Du et al. (2016).
- `machine-learning-scoring-functions.md`: Added citations for Huang et al. (2010) and van Pinxteren & Jespers (2025).
- `mm-pbsa-and-mm-gbsa.md`: Added citation for Bruce et al. (2018).
- `molecular-dynamics-simulation.md`: Added citations for Salmaso & Moro (2018), Bernardi et al. (2015), Genheden & Ryde (2015), and Poli et al. (2020).
- `protein-ligand-docking.md`: Added citations for Eberhardt et al. (2021) and Salmaso & Moro (2018).
- `protein-ligand-interactions.md`: Added citations for Bruce et al. (2018) and Eberhardt et al. (2021).
- `receptor-preparation.md`: Added citation for Hospital et al. (2015).
- `scoring-functions.md`: Added citations for Yang et al. (2022) and Du et al. (2016).
- `virtual-screening.md`: Added citations for Huang et al. (2010), Eberhardt et al. (2021), and Li et al. (2025).

**Contradiction and outdated claim updates (Findings #4, #5, #6, #7)**

- *Finding #4 (CELPP vs Self-Docking)*: Differentiated retrospective redocking (70–85% pose reproduction in pre-formed cavities) from prospective blind scoring across `scoring-functions.md`, `protein-ligand-docking.md`, `benchmarking-datasets.md`, and `huang-2010-scoring-functions.md`, highlighting the CELPP finding that conformational sampling succeeded in 80.4% of targets while scoring correctly ranked the native pose first in only 46.3%.
- *Finding #5 (PoseBusters Validity)*: Explicitly defined PoseBusters (PB-Valid) metrics (clash ratio < 0.075, bond length/angle tolerances, planarity, tetrahedral chirality) across `deep-learning-docking.md`, `machine-learning-scoring-functions.md`, and `protein-ligand-docking.md`, detailing KarmaDock's 0.00% PB-valid rate across Astex/DockGen and DiffDock's receptor steric clash frequencies requiring downstream force-field post-relaxation.
- *Finding #6 (AutoDock Vina MW Bias and Metalloenzyme Handling)*: Contextualized Vina's empirical scoring function in `autodock-vina.md` and `vina-docs-basic-docking.md` using Xu et al.'s benchmark, detailing the molecular weight bias penalizing fragments, generation of 311 irregular torsions, and the requirement for explicit partial charges or AutoDock4Zn in metalloenzymes like COMT.
- *Finding #7 (Homology and AlphaFold Viability in FEP)*: Contextualized the older assumption that FEP strictly requires <2.5 Å crystal structures across `force-fields.md`, `binding-free-energy.md`, and `alchemical-free-energy-methods.md`, documenting Muegge & Hu (2023) and van Pinxteren & Jespers (2025) evidence that homology models (30–50% identity) and AlphaFold holo models achieve crystal-level accuracy ($\text{RMSE} \sim 1.0\text{--}1.5\text{ kcal/mol}$) following initial equilibrium MD relaxation.

**Navigation and index integration**

- Updated `wiki/index.md` registering `experimental-binding-assays.md`, `hydrogen-mass-repartitioning.md`, and `benchmark-drug-targets.md` with one-line descriptions.
- Re-anchored reading paths for Intermediate (Workflows & Failures), Advanced (Free Energy, Kinetics & DL), and MD Specialist to incorporate the new concept pages.
- Total wiki pages: **53 markdown pages** (28 concept pages, 23 source summaries, `index.md`, `log.md`).


---

## 2026-11-01 — Roadmap completion: scientific hedging, citation audit, secondary link graph reinforcement, and four new concept pages

**Context and scope**

Executed the complete 5-step prioritized remediation roadmap following the comprehensive wiki audit to achieve total structural, physical, and topological integrity across the entire knowledge base according to `AGENTS.md` guidelines.

**Step 1: Scientific hedging and contradiction updates (Findings #1, #2, #3, #4, #5)**

- *Resolution threshold hedging (Finding #1)*: Updated `molecular-dynamics-simulation.md`, `structure-based-drug-design.md`, and `decherchi-2020-binding-thermodynamics-kinetics.md` to reconcile the classical rule of thumb requiring <2.5 Å crystal structures with modern evidence from Muegge & Hu (2023) and van Pinxteren & Jespers (2025), noting that homology models (30–50% identity) and AlphaFold/Chai-1/Boltz-1 models routinely achieve crystal-grade accuracy after equilibrium MD relaxation.
- *Continuum electrostatics solvation nuances (Finding #2)*: Updated `binding-free-energy.md` to reconcile continuum electrostatics ranking performance citing Genheden & Ryde (2015), documenting that while PB represents more rigorous classical continuum physics, GB often equals or outperforms PB in empirical ranking benchmarks due to systematic cancellation of dielectric boundary discretization errors.
- *Water molecules handling and hydrated docking (Finding #3)*: Updated `receptor-preparation.md` to clarify the operational distinction between stripping bulk, disordered crystallographic water molecules during initial receptor preparation and preserving structurally conserved, bridging waters with negative hydration free energies (via WaterMap/3D-RISM) for hydrated docking and alchemical FEP.
- *AI co-folding model citations (Finding #4)*: Updated `binding-site-detection.md`, `protein-ligand-docking.md`, `receptor-preparation.md`, and `virtual-screening.md` to explicitly cite Morehead et al. (2026) for modern all-atom biomolecular co-folding architectures (AlphaFold 3, Chai-1, Boltz-1).
- *Virtual screening engine throughput contextualization (Finding #5)*: Updated `virtual-screening.md` and `autodock-vina.md` to contextualize AutoDock Vina's throughput (100-fold faster than AutoDock 4, but approximately 23-fold slower than rigid rigid-fragment grids like DOCK 3.7 on DUD-E benchmark screening campaigns) citing Xu et al. (2022).

**Step 2: Universal inline citation anchoring for quantitative metrics (Finding #11)**

- Added explicit inline `(source: ...)` citations to all 18 quantitative bullet points identified in the audit across `machine-learning-scoring-functions.md`, `virtual-screening.md`, `binding-site-detection.md`, `binding-free-energy.md`, `binding-thermodynamics.md`, `binding-kinetics.md`, `ligand-preparation-and-protonation.md`, `receptor-preparation.md`, and `benchmarking-datasets.md`.
- Automated verification confirms zero unanchored quantitative claims across all 32 concept pages.

**Step 3: Secondary inbound link graph reinforcement (Finding #10)**

- Reinforced network topology by adding secondary inbound wiki-links to all 14 previously single-inbound source summary pages (`ma-2024-docking-strategies`, `huang-2010-scoring-functions`, `yang-2022-ml-docking`, `durrant-2011-md-drug-discovery`, `hollingsworth-2018-md-for-all`, `du-2016-protein-ligand-interactions`, `vina-docs-basic-docking`, `salmaso-2018-bridging-docking-md`, `xu-2022-docking-failures-analysis`, `genheden-2015-mm-pbsa-gbsa`, `poli-2020-mm-pbsa-virtual-screening`, `bernardi-2015-enhanced-sampling-md`, `lazim-2020-md-enhanced-sampling-advances`, `morehead-2026-evaluating-dl-docking`).
- Every source summary page now possesses $\ge 2$ inbound links from relevant concept pages.

**Step 4: Creation and cross-linking of four new concept pages**

1. `metadynamics.md` — History-dependent repulsive bias potential deposition along collective variables ($V(s, t) = \sum_{t'} W \exp\left(-\sum_{i} \frac{(s_i - s_i(t'))^2}{2\sigma_i^2}\right)$), Well-Tempered Metadynamics bias factor rescaling ($\Delta T / (T + \Delta T)$), Infrequent Metadynamics (InMetaD) unbinding kinetics ($k_{\text{off}}$), and transition acceleration factors ($\alpha$).
2. `umbrella-sampling.md` — Stratified equilibrium sampling along reaction coordinates using harmonic bias potentials ($V_i(\xi) = \frac{1}{2} k_i (\xi - \xi_{0,i})^2$), overlap criteria, WHAM and MBAR unbiasing, PMF integration, and standard-state volume corrections ($V_0 = 1661\text{ \AA}^3$).
3. `markov-state-models.md` — Master-equation kinetic modeling from distributed MD ensembles ($\mathbf{p}(t + \tau) = \mathbf{p}(t) \mathbf{T}(\tau)$), TICA slow coordinate projection, $k$-means microstate partitioning, Chapman–Kolmogorov test implied timescales validation, and Transition Path Theory (TPT) reactive flux.
4. `chemical-screening-artifacts.md` — Systematic taxonomy of promiscuous false-positive mechanisms in virtual and experimental screening: PAINS reactive chemotypes (rhodanines, quinones, curcumin), colloidal aggregator self-assembly ($100\text{--}1000\text{ nm}$ particles), Triton X-100 detergent counter-screens ($0.01\%\text{ v/v}$), and orthogonal biophysical validation (ITC, SPR, 2D HSQC NMR).
- Fully cross-linked the 4 new pages into `enhanced-sampling.md`, `binding-kinetics.md`, `binding-free-energy.md`, `benchmark-drug-targets.md`, `cryptic-and-allosteric-pockets.md`, `virtual-screening.md`, `docking-failure-modes.md`, and `experimental-binding-assays.md`.

**Step 5: Index and navigation updates**

- Updated `wiki/index.md` registering the 4 new concept pages with one-line descriptions.
- Re-anchored `## Reading paths` across Intermediate (Workflows & Failures), Advanced (Free Energy, Kinetics & DL), and MD Specialist to integrate `[[chemical-screening-artifacts]]`, `[[markov-state-models]]`, `[[umbrella-sampling]]`, and `[[metadynamics]]`.
- Total wiki pages: **57 markdown files** (32 concept pages, 23 source summaries, `index.md`, `log.md`).
- Wiki health: 0 broken links, 0 orphaned pages, 100% adherence to AGENTS.md template and citation rules.


---

## 2026-11-02 — Systematic ingestion of 7 new molecular simulation and benchmark sources

**Context and scope**

Scanned the source repository `row/` and identified seven unlogged foundational files covering the GROMACS protein–ligand simulation protocol, LiveCoMS simulation and affinity benchmark best practices, Mixed-Solvent MD (MixMD), direct chemical perception and Open Force Field (Sage), continuous neural network force fields (Garnet), and MD trajectory time-series descriptors. Ingested all seven sources systematically according to `AGENTS.md` guidelines, creating seven dedicated source summary pages, four new concept/protocol pages, cross-linking existing pages, and updating navigation.

**Sources read and ingested (7)**

1. `gromacs_protein_ligand_complex_tutorial_lemkul.md` — Justin A. Lemkul's complete step-by-step tutorial on setting up, parameterizing, equilibrating, simulating, and analyzing a protein–ligand complex (T4 lysozyme L99A/M102Q bound to 2-propylphenol, JZ4) using CHARMM36 and CGenFF in GROMACS.
2. `mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf` — Lemkul's peer-reviewed LiveCoMS article presenting a standardized suite of seven GROMACS tutorials (lysozyme, membrane proteins, umbrella sampling, biphasic systems, ligands, alchemical free energy of solvation, and virtual sites).
3. `Best Practices for Foundations in Molecular Simulations [Article.pdf` — Braun et al. (LiveCoMS, 2019): comprehensive foundational best practices, statistical mechanical ensembles, integrator and timestep criteria, thermostat and barostat evaluations (flying ice cube artifact, Bussi v-rescale, Nosé–Hoover), PBC conventions, and setup checklists.
4. `Best practices for constructing, preparing, and evaluating protein-ligand binding affinity.pdf` — Hahn et al. (LiveCoMS, 2022): guidelines for binding affinity benchmarks (PLBenchmarks and arsenic), biophysical assay curation (ITC/SPR vs. enzymatic IC50), experimental noise floors ($\sim 0.3\text{--}0.5\text{ kcal/mol}$), crystal quality filtering ($R_{\text{free}} < 0.25$, RSCC), and closed-loop perturbation networks.
5. `The Open Force Field Initiative Open Software and Open Science.pdf` — Wang et al. (J. Phys. Chem. B, 2024): the Open Force Field Initiative, direct chemical perception via SMIRNOFF, eliminating indirect atom-type explosion, automated parameterization via ForceBalance and QCArchive, and the evolution of open-source force fields from Parsley to Sage 2.2.0.
6. `Training a force field for proteins and small molecules from.pdf` — Blanco-González et al. (MRC LMB, 2024): Garnet, an end-to-end graph neural network parameterizing all molecular mechanics terms from scratch for both proteins and small molecules using continuous atom typing and a double exponential nonbonded potential.
7. `Assessment of molecular dynamics time series.pdf` — Poziemski et al. (Digital Discovery / RSC, 2026): large-scale benchmark on 2,752 complexes evaluating machine learning models on MD time-series descriptors (`ts_descriptors`), demonstrating that multi-stage feature trimming isolates contact persistence and energy fluctuations that overcome static crystal limitations on target and scaffold generalization splits.

**Pages created — source summaries (7)**

- `lemkul-gromacs-protein-ligand-tutorial.md`
- `lemkul-2018-gromacs-tutorials.md`
- `braun-2019-foundations-molecular-simulation.md`
- `hahn-2022-binding-affinity-benchmarks.md`
- `wang-2024-open-force-field-initiative.md`
- `blanco-gonzalez-2024-garnet-force-field.md`
- `poziemski-2026-md-time-series-descriptors.md`

**Pages created — concepts and protocols (4)**

1. `gromacs-protein-ligand-protocol.md` — Synthesizes the end-to-end GROMACS protocol: coordinate cleaning, protein topology generation with `pdb2gmx`, small-molecule ligand parameterization via CGenFF (`.mol2` bond sorting, penalty score validation rules $<10, 10\text{--}50, >50$), rhombic dodecahedral box definition (saving $\sim 23\%$ solvent volume), ion neutralization, CHARMM36 force-switch nonbonded settings (`rvdw-switch = 1.0`, `rvdw = 1.2`), ligand heavy-atom position restraints, merged temperature coupling groups (`tc-grps = Protein_JZ4 Water_and_ions`), LINCS constraints, production MD, PBC re-wrapping via `gmx trjconv -pbc mol -ur compact`, hydrogen-bond geometry, heavy-atom ligand RMSD, and short-range interaction energy decomposition via `mdrun -rerun`.
2. `mixed-solvent-molecular-dynamics.md` — Complete technical guide to MixMD / cosolvent MD (MixMD, SILCS, SWISH, CrypticScout): organic probe selection spanning pharmacophoric features (isopropanol, benzene, phenol, acetonitrile, acetate, methylammonium), concentration thresholds (5% v/v to avoid phase separation/denaturation), multiple replicate convergence, 3D Grid Free Energy (GFE) calculation via Boltzmann inversion ($\Delta G_{\text{probe}}(\mathbf{r}) = -k_B T \ln [\rho(\mathbf{r}) / \rho_{\text{bulk}}]$), and translation into Solvent-Site-Based Docking (SSBD) and dynamic pharmacophores.
3. `open-force-field-and-ligand-parameterization.md` — Detailed comparative analysis of small-molecule parameterization: indirect atom typing and atom-type explosion in GAFF/CGenFF/OPLS vs. direct chemical perception via SMIRNOFF SMARTS/SMIRKS graph patterns in OpenFF Sage 2.2.0; continuous atom typing and end-to-end neural force fields in Garnet utilizing the double exponential potential ($U_{\text{dexp}}(r)$) to ensure stable optimization gradients on GPUs.
4. `simulation-and-benchmark-best-practices.md` — Consolidated best-practices manual merging LiveCoMS simulation guidelines (integrators, timesteps $\Delta t \le \tau_{\text{fastest}}/10$, Bussi/Nosé–Hoover thermostats, Parrinello–Rahman barostats, PBC minimum image convention $L \ge 2 r_{\text{cut}}$, replicate sampling) with Hahn et al.'s standards for binding affinity benchmarks (biophysical ITC/SPR assay curation, avoiding raw IC50 conversions, structural quality metrics, closed thermodynamic cycles, and avoiding misleading $R^2$ inflation).

**Existing pages updated and cross-linked (9)**

- `molecular-dynamics-simulation.md`: Integrated LiveCoMS simulation checklists, Bussi vs. Berendsen thermostat artifacts, rhombic dodecahedral cell volume savings, Lemkul GROMACS protocol, MixMD cryptic site discovery, and MD time-series descriptors.
- `force-fields.md`: Integrated OpenFF direct chemical perception (SMIRNOFF, Sage 2.2.0), Garnet continuous GNN force field and double exponential potential, and CGenFF penalty validation rules.
- `ligand-preparation-and-protonation.md`: Added dedicated small-molecule parameterization workflow for MD (CGenFF, GAFF/ACPYPE, OpenFF Sage, Garnet) and benchmark assay pH state curation.
- `receptor-preparation.md`: Integrated Hahn et al. crystallographic benchmark criteria (resolution $\le 2.0\text{--}2.5\text{ \AA}$, $R_{\text{free}} < 0.25$, RSCC), Lemkul's `pdb2gmx` setup, and cosolvent probe mapping.
- `binding-site-detection.md`: Added explicit cosolvent MD (MixMD / SILCS) probe mapping, 3D GFE calculations, and solvent-site-based docking.
- `cryptic-and-allosteric-pockets.md`: Added explicit cosolvent MD probe induction mechanisms and GFE hotspot mapping.
- `benchmarking-datasets.md`: Added PLBenchmarks and arsenic, the 2,752-complex MD time-series dataset, FreeSolv, and SPICE quantum mechanical datasets.
- `binding-free-energy.md`: Integrated Hahn et al. benchmark design criteria (closed perturbation cycles, cycle closure hysteresis, ITC/SPR curation, experimental error floor) and OpenFF Sage / Garnet RBFE benchmark performance.
- `machine-learning-scoring-functions.md`: Added dedicated section on MD time-series descriptors (`ts_descriptors`), multi-stage feature trimming, and Garnet differentiable force-field learning.

**Index and graph verification**

- Registered all 11 new pages (4 concepts + 7 source summaries) in `wiki/index.md` and integrated them into the **Advanced** and **MD Specialist** reading paths.
- Total wiki pages: **68 markdown files** (36 concept pages, 30 source summaries, `index.md`, and `log.md`).
- Graph integrity: **0 broken links**, **0 orphaned pages**, every new page has $\ge 3$ inbound links from relevant content pages.
