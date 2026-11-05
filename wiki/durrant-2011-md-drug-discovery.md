# Durrant & McCammon (2011) — MD Simulations and Drug Discovery

**Summary**: A review arguing that molecular dynamics fills the gap between static structures and the dynamic reality of drug binding, with applications in cryptic-site discovery, the relaxed complex scheme, and alchemical free-energy calculations.
**Sources**: Molecular dynamics simulations and drug discovery.md
**Last updated**: 2026-11-01

---

## Citation and scope

Jacob D. Durrant and J. Andrew McCammon, "Molecular dynamics simulations and drug discovery," *BMC Biology* 9:71 (2011), doi:10.1186/1741-7007-9-71. PMC3203851 (source: Molecular dynamics simulations and drug discovery.md).

The review is written for a broad biology audience. It explains MD in physical terms, reviews its limitations, and then surveys three drug-discovery applications: cryptic/allosteric site discovery, the relaxed complex scheme, and alchemical free-energy calculations.

## The argument

A single crystal structure is like a photograph of a runner: it says little about the stride. Proteins are in constant motion, and both molecular recognition and drug binding are dynamic processes. MD provides the movie rather than the snapshot (source: Molecular dynamics simulations and drug discovery.md).

The review is careful not to oversell the method. Its two central limitations are force-field accuracy and computational cost, and it notes that alchemical free-energy calculations are uniquely sensitive to inadequate sampling (source: Molecular dynamics simulations and drug discovery.md).

## MD in plain terms

MD was first developed in the late 1970s. It replaces quantum-mechanical motion with Newtonian physics: bonds and angles are springs, dihedrals are sinusoidal functions, van der Waals interactions use the Lennard-Jones potential, and electrostatics use Coulomb's law. The parameters collectively form a force field, fitted to quantum-mechanical calculations and experimental data (source: Molecular dynamics simulations and drug discovery.md).

Positions are updated by Newton's laws, the time is advanced by 1–2 fs, and the process repeats millions of times. Simulations are usually run on clusters or supercomputers with MPI parallelism (source: Molecular dynamics simulations and drug discovery.md).

Common force fields are AMBER, CHARMM, and GROMOS; they differ mainly in parameterization but generally give similar results (source: Molecular dynamics simulations and drug discovery.md).

## Limitations identified in 2011

- **Force fields need refinement.** Classical MD cannot break or form bonds, ignores electronic polarization, and handles transition metals poorly. QM/MM addresses the first by treating a small region quantum mechanically (source: Molecular dynamics simulations and drug discovery.md).
- **Computational cost.** A 1 µs simulation of about 25,000 atoms on 24 processors took several months. Most simulations were nanoseconds, and many processes of interest occur on much longer timescales (source: Molecular dynamics simulations and drug discovery.md).
- **Sampling.** Insufficient sampling can miss relevant conformations; the results are then incomplete rather than necessarily wrong (source: Molecular dynamics simulations and drug discovery.md).

Hardware responses discussed: accelerated MD (aMD) lowers barriers, GPUs speed simulations by about an order of magnitude, and Anton was designed specifically for MD, reaching microseconds per day and simulations longer than 1 ms that captured folding, unfolding, and drug binding (source: Molecular dynamics simulations and drug discovery.md).

## Application 1: cryptic and allosteric sites

The review defines cryptic sites as pockets not obvious in available structures and describes MD as an excellent tool for finding them (source: Molecular dynamics simulations and drug discovery.md).

The flagship example is HIV integrase. In 2004, Schames and coworkers simulated the protein and found a trench absent from all crystal structures. Crystallography later showed that known inhibitors bind there. The work led to new studies at Merck, further development, and ultimately **raltegravir**, the first FDA-approved HIV integrase inhibitor (source: Molecular dynamics simulations and drug discovery.md).

A second example is allosteric site discovery on β1 and β2 adrenergic receptors. Ivetac and McCammon extracted multiple conformations and flooded the surface with organic probes using FTMAP; regions where probes congregated across structures were nominated as allosteric sites. Five sites were found, some absent from crystal structures (source: Molecular dynamics simulations and drug discovery.md).

## Application 2: the relaxed complex scheme

Traditional docking uses one static receptor, which can discard true ligands that bind a different conformation. The relaxed complex scheme (RCS) docks each ligand into multiple MD-derived conformations, giving a spectrum of scores per ligand; ligands are ranked by a statistic such as the mean score (source: Molecular dynamics simulations and drug discovery.md).

RCS has identified inhibitors of FKBP, HIV integrase, *Trypanosoma brucei* RNA editing ligase 1, *T. brucei* GalE, *T. brucei* FPPS, and *Mycobacterium tuberculosis* dTDP-6-deoxy-L-lyxo-4-hexulose synthase. In two projects the inhibitors were effective against the whole parasite (source: Molecular dynamics simulations and drug discovery.md).

Weaknesses: RCS inherits force-field and sampling errors from MD, and it relies on docking scoring functions that sacrifice accuracy for speed, treating conformational entropy and solvation superficially (source: Molecular dynamics simulations and drug discovery.md).

## Application 3: alchemical free energy

Free energy is a state function, so the path from unbound to bound does not matter. This licenses alchemical transformations in which the ligand's electrostatic and van der Waals interactions are gradually turned off — the ligand "disappears" — in the binding site and in solution (source: Molecular dynamics simulations and drug discovery.md).

The thermodynamic cycle gives:

$$\Delta G_{bind} = \Delta G_{water} - \Delta G_{protein}$$

where ΔG_protein is the free energy of annihilating the bound ligand and ΔG_water is the free energy of annihilating the solvated ligand. The "ghost" ligand binding step is zero because the ghost does not interact (source: Molecular dynamics simulations and drug discovery.md).

Relative binding free energies are computed by mutating one part of a ligand into another — for example, converting a carbon to an oxygen — rather than annihilating the whole molecule (source: Molecular dynamics simulations and drug discovery.md).

The review notes a cycle of enthusiasm and disappointment: early FEP results in the 1980s and 1990s agreed remarkably with experiment, but subsequent predictions fell short. Algorithmic and engineering advances renewed interest, with accurate predictions for the src SH2 domain, a mutant T4 lysozyme, FKBP12, HIV reverse transcriptase, trypsin, a bacterial ribosome, and estrogen receptor α (source: Molecular dynamics simulations and drug discovery.md).

The main caveat is sampling. Alchemical calculations are far more sensitive to incomplete sampling than RCS screens, because unsampled conformations do not contribute to the calculated free energy. This is why they have not been widely adopted by industry despite great interest (source: Molecular dynamics simulations and drug discovery.md).

## Key numbers

| Quantity | Value |
|---|---|
| 1 µs, ~25,000 atoms, 24 CPUs | several months |
| Typical timestep | 1–2 fs |
| Anton | microseconds per day; >1 ms simulations |
| RCS inhibitors validated | FKBP, HIV integrase, Tb ligase 1, Tb GalE, Tb FPPS, Mtb enzyme |
| Allosteric sites found on β-AR | 5 |

(source: Molecular dynamics simulations and drug discovery.md)

## Strengths

- Exceptionally clear physical explanations of MD, force fields, and alchemical free energy.
- Uses concrete drug-discovery success stories, especially raltegravir.
- Honest about limitations and the history of over-optimism.
- Introduces the relaxed complex scheme in accessible terms.

## Limitations

- Published in 2011; hardware, force fields, and enhanced sampling have advanced substantially. Its cost figures are historical, not current.
- Focuses on the authors' interests (RCS, alchemical methods) and does not review enhanced sampling in depth.
- The raltegravir story is presented as a success; later work has refined the details of integrase inhibition.
- No quantitative benchmark of MD accuracy.

## Why it matters for this wiki

This is the source for [[cryptic-and-allosteric-pockets]], [[molecular-dynamics-simulation]], [[binding-free-energy]], and the relaxed complex scheme discussion in [[virtual-screening]]. It supplies the raltegravir case study and the conceptual basis for alchemical free energy.

## Related pages
- [[molecular-dynamics-simulation]]
- [[cryptic-and-allosteric-pockets]]
- [[binding-free-energy]]
- [[enhanced-sampling]]
- [[virtual-screening]]
- [[hollingsworth-2018-md-for-all]]
