# Gimeno et al. (2019) — The Light and Dark Sides of Virtual Screening

**Summary**: A critical review detailing the methodologies, advantages, and deceptive pitfalls of ligand-based and structure-based virtual screening in drug discovery.
**Sources**: The Light and Dark Sides of Virtual Screening.pdf
**Last updated**: 2026-11-01

---

## Citation and scope

Aleix Gimeno, María José Ojeda-Montes, Sarah Tomás-Hernández, Adrià Cereto-Massagué, Raúl Beltrán-Debón, Miquel Mulero, Gerard Pujadas, and Santiago Garcia-Vallvé, "The Light and Dark Sides of Virtual Screening: What Is There to Know?" *International Journal of Molecular Sciences* 20(6):1375 (2019), doi:10.3390/ijms20061375 (source: The Light and Dark Sides of Virtual Screening.pdf).

The review examines the end-to-end virtual screening (VS) pipeline, juxtaposing computational strengths ("the light side") against prevalent experimental and algorithmic failure modes ("the dark side") that cause false positives and wasted resources in hit discovery.

## The light side: screening efficiency and methodologies

Virtual screening processes tens of thousands to millions of virtual or purchasable molecules in hours to days, reducing wet-lab screening costs and expanding accessible chemical space beyond physical High-Throughput Screening (HTS) inventories (source: The Light and Dark Sides of Virtual Screening.pdf).

The review systematically compares two foundational computational methodologies:

### 1. Ligand-based virtual screening (LBVS)
LBVS operates on the similar property principle formulated by Johnson and Maggiora: structurally similar compounds display similar biological effects (source: The Light and Dark Sides of Virtual Screening.pdf). Key approaches include:
- **Molecular fingerprints and similarity**: 2D topological path fingerprints, circular fingerprints (ECFP, FCFP), and pharmacophore fingerprints evaluated via Tanimoto, Dice, or Cosine metrics (source: The Light and Dark Sides of Virtual Screening.pdf).
- **Machine learning classification**: Random Forest (RF), Support Vector Machines (SVM), $k$-Nearest Neighbors ($k$-NN), and Artificial Neural Networks (ANN) trained on known actives and inactives (source: The Light and Dark Sides of Virtual Screening.pdf).
- **Shape similarity and 3D pharmacophores**: Alignment of ligand steric shapes and electrostatic features using tools such as ROCS or SHAEP to identify scaffold-hopping analogs (source: The Light and Dark Sides of Virtual Screening.pdf).

### 2. Receptor-based virtual screening (RBVS / SBVS)
RBVS docks candidate molecules into a three-dimensional receptor pocket derived from X-ray crystallography, NMR, cryo-EM, or homology modeling (source: The Light and Dark Sides of Virtual Screening.pdf). The computational process entails:
- Protein preparation (protonation, tautomeric states, loop repair, and handling of structural water molecules and metal ions);
- Cavity and binding site detection;
- Conformational sampling via incremental construction or stochastic search (Monte Carlo, genetic algorithms);
- Affinity estimation using force-field, empirical, or knowledge-based scoring functions (source: The Light and Dark Sides of Virtual Screening.pdf).

## The dark side: false positives and artifacts

The central warning of the review is that top-scoring computational hits frequently prove inactive or artifactual in experimental validation due to specific, well-characterized mechanisms (source: The Light and Dark Sides of Virtual Screening.pdf):

### 1. Pan-assay interference compounds (PAINS)
PAINS represent chemical sub-structures (e.g., rhodanines, hydroxyphenyl hydrazones, curcuminoids, toxoflavins) that trigger false-positive signals across multiple unrelated biological assays (source: The Light and Dark Sides of Virtual Screening.pdf). Interference mechanisms include:
- Covalent modification of protein residues;
- Redox activity generating hydrogen peroxide ($H_2O_2$);
- Direct fluorescence emission or fluorescence quenching;
- Metal chelation disrupting metalloenzymes (source: The Light and Dark Sides of Virtual Screening.pdf).
The authors emphasize filtering candidate libraries against PAINS patterns (e.g., via PAINS-Remover) prior to or immediately following docking.

### 2. Colloidal aggregators
Many hydrophobic, promiscuous small molecules aggregate at micromolar concentrations in aqueous buffer to form colloidal particles 100 to 1,000 nm in diameter (source: The Light and Dark Sides of Virtual Screening.pdf). These aggregates non-specifically sequester, unfold, or inhibit enzymes without specific active-site binding (source: The Light and Dark Sides of Virtual Screening.pdf).
- **Diagnostic criteria**: Identification via computational aggregators tools (Aggregator Advisor);
- **Experimental verification**: Inhibition sensitivity to non-ionic detergents (e.g., 0.01%–0.1% Triton X-100) or enzyme concentration changes (source: The Light and Dark Sides of Virtual Screening.pdf).

### 3. Decoy bias in benchmark validation
The authors caution against over-optimism when evaluating VS protocols on synthetic benchmarks like DUD and DUD-E. While DUD-E provides 50 property-matched decoys per active to control for charge, molecular weight, and $\log P$, algorithms often exploit residual physical biases or unrepresentative active clustering rather than true molecular recognition signatures (source: The Light and Dark Sides of Virtual Screening.pdf).

## Practical recommendations for hit selection

To maximize genuine hit rates and minimize experimental attrition, the authors outline an integrated hit selection protocol (source: The Light and Dark Sides of Virtual Screening.pdf):
1. **Combine ranking with visual inspection**: Eliminate docking poses that place hydrophobic moieties into bulk solvent or leave buried polar residues without counter-hydrogen bonds.
2. **Interaction fingerprints (SIFt)**: Filter docking poses based on the presence of essential native interactions with catalytic residues.
3. **Scaffold clustering**: Cluster top-scoring poses to select chemically diverse chemotypes rather than many analogs of a single high-scoring core.
4. **Physicochemical and ADMET filters**: Enforce Lipinski's Rule of 5 and Veber rules early to eliminate intractable compounds.
5. **Orthogonal biophysical assays**: Validate hits using surface plasmon resonance (SPR), isothermal titration calorimetry (ITC), or NMR to confirm stoichiometric, non-aggregated binding (source: The Light and Dark Sides of Virtual Screening.pdf).

## Related pages
- [[virtual-screening]]
- [[protein-ligand-docking]]
- [[scoring-functions]]
- [[ligand-preparation-and-protonation]]
- [[receptor-preparation]]
- [[benchmarking-datasets]]

