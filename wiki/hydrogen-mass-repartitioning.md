# Hydrogen Mass Repartitioning

**Summary**: A simulation method that redistributes mass from heavy atoms to covalently bound hydrogens, slowing high-frequency vibrational motions and enabling stable 4.0 fs integration timesteps in molecular dynamics.
**Sources**: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Molecular dynamics simulation for all.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md
**Last updated**: 2026-11-01

---

## Physical motivation and numerical integration limits

In classical [[molecular-dynamics-simulation]], the maximum allowable integration timestep ($\Delta t$) is governed by the highest-frequency vibrational motions in the simulated system (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Molecular dynamics simulation for all.md). According to the Verlet algorithm and symplectic integrators, numerical stability requires that the integration timestep does not exceed approximately one-tenth of the shortest vibrational period (source: Molecular dynamics simulation for all.md).

For biological macromolecules, the fastest degrees of freedom are the covalent bond-stretching vibrations involving hydrogen atoms (such as $\text{C}-\text{H}$, $\text{N}-\text{H}$, and $\text{O}-\text{H}$ bonds), which possess vibrational frequencies on the order of $3000\text{--}3500\text{ cm}^{-1}$ and oscillation periods of approximately $10\text{ fs}$ (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Molecular dynamics simulation for all.md):
- Without bond constraints, simulations are restricted to a timestep of $\Delta t = 1.0\text{ fs}$ to avoid numerical instability and energy drift (source: Molecular dynamics simulation for all.md).
- Applying holonomic constraints (such as SHAKE, RATTLE, or LINCS) to all covalent bonds involving hydrogen removes bond-stretching frequencies, extending the stable integration timestep to $\Delta t = 2.0\text{ fs}$ (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Molecular dynamics simulation for all.md).
- However, extending beyond $2.0\text{ fs}$ in standard topologies causes high-frequency angle-bending modes involving light hydrogen atoms (oscillation period $\sim 20\text{ fs}$) to resonate with the integrator, leading to rapid kinetic energy accumulation and system instability (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

## The hydrogen mass repartitioning (HMR) protocol

Hydrogen Mass Repartitioning (HMR) overcomes the $2.0\text{ fs}$ threshold by altering the inertial properties of the solute without altering its potential energy surface (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

### 1. Mass transfer formulation
HMR reallocates atomic mass from heavy parent atoms (carbon, nitrogen, oxygen, sulfur) to their covalently bonded hydrogen atoms (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf):
- The mass of each hydrogen atom ($m_\text{H} \approx 1.008\text{ amu}$) is increased by a scaling factor, typically to **$\sim 3.0\text{ amu}$** (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).
- An equivalent total mass ($\Delta m = \sum (m_\text{H, new} - m_\text{H, standard})$) is deducted from the bonded heavy parent atom ($m_\text{heavy, new} = m_\text{heavy, standard} - \Delta m$) (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).
- This formulation ensures that the **total mass of the residue, molecule, and system remains strictly invariant**, and the center of mass remains conserved (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

### 2. Frequency lowering
Because the harmonic vibrational frequency of an angle bend or bond oscillation scales inversely with the square root of the reduced mass ($\omega = \sqrt{k/\mu}$), tripling the mass of hydrogen decreases the high-frequency bending modes below the resonance limit (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). Coupled with holonomic constraints on all bonds (e.g., LINCS or SHAKE), this modification enables numerically stable integration at:
$$\Delta t = 4.0\text{ fs}$$
doubling the simulation throughput per unit of wall-clock GPU time (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

## Thermodynamic and kinetic validation

Because the potential energy function $V(\mathbf{r})$ is completely independent of atomic masses, the configuration-space partition function $Z_N$ and all equilibrium thermodynamic properties in the canonical ($NVT$) and isothermal-isobaric ($NPT$) ensembles are theoretically invariant under mass alteration (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Empirical validation benchmarks across diverse biological systems confirm this theoretical invariance:
- **Peptides and globular proteins**: Hopkins et al. evaluated HMR across benchmark systems ranging from tripeptides to the 129-residue Hen Egg White Lysozyme (HEWL) (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). The authors confirmed that conformational equilibrium ensembles, backbone root-mean-square fluctuations ([[molecular-dynamics-simulation]]), and free energy landscapes were indistinguishable between $2.0\text{ fs}$ standard mass and $4.0\text{ fs}$ HMR trajectories (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).
- **Lipid bilayers and membrane proteins**: Balusek et al. and Wang & Gumbart assessed HMR validity in complex lipid membrane simulations (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). Equilibrium structural parameters—including area per lipid, electron density profiles, and deuterium order parameters ($S_\text{CD}$)—exhibited complete agreement with standard simulations (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf). Subtle kinetic variations, such as modest reductions in lateral lipid diffusion coefficients, were noted due to inertial dampening, but internal conformational dynamics remained preserved (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

## Integration with simulation workflows

### 1. Water models
Standard rigid 3-site and 4-site water models (e.g., TIP3P, SPC/E, TIP4P-Ew) utilize analytical constraint algorithms like SETTLE to hold all O-H bonds and H-O-H angles rigid (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Molecular dynamics simulation for all.md). Consequently, water molecules do not suffer from flexible angle-bending resonance at $4.0\text{ fs}$ and generally retain standard masses ($m_\text{H} \approx 1.008\text{ amu}$, $m_\text{O} \approx 16.00\text{ amu}$), although solute macromolecules (proteins, nucleic acids, lipids, and docked ligands) undergo HMR repartitioning (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

### 2. Free energy calculations and enhanced sampling
- In alchemical free energy calculations ([[alchemical-free-energy-methods]]), HMR is widely adopted in both relative binding free energy (RBFE) and absolute binding free energy (ABFE) workflows to accelerate intermediate alchemical $\lambda$-window sampling (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- In enhanced sampling ([[enhanced-sampling]]), HMR is fully compatible with metadynamics, accelerated MD (aMD), and adaptive sampling, effectively doubling phase-space exploration within fixed compute budgets (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- In [[binding-kinetics]] estimation, when extracting absolute unbinding rate constants ($k_\text{off}$) or residence times, unadjusted dynamic rates may show minor inertial friction shifts; authors recommend reporting kinetic constants with explicit notation of mass repartitioning parameters (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Related pages
- [[molecular-dynamics-simulation]]
- [[force-fields]]
- [[enhanced-sampling]]
- [[alchemical-free-energy-methods]]
- [[binding-free-energy]]
- [[binding-kinetics]]
- [[lazim-2020-md-enhanced-sampling-advances]]
