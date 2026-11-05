# Braun et al. (2019) — Foundations in Molecular Simulations

**Summary**: The LiveCoMS foundational best practices guide establishing essential physical, mathematical, and algorithmic guidelines for molecular dynamics simulations, featuring a comprehensive setup checklist and a comparative audit of thermostats and barostats.
**Sources**: `Best Practices for Foundations in Molecular Simulations [Article.pdf`
**Last updated**: 2026-11-02

---

## Foundational principles and simulation design

Braun et al. provide a systematic starting point for molecular simulations, emphasizing that accurate biomolecular modeling requires deep alignment between the physical question, the thermodynamic ensemble, and the numerical algorithm (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). Classical Molecular Dynamics (MD) propagates particles by integrating Newton's equations of motion across a potential energy surface defined by a molecular mechanics (MM) force field (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). Because classical force fields hold bonding topology invariant (barring reactive formulations), chemistry and protonation states remain fixed throughout the trajectory (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

The guide emphasizes the principle of "counting the cost": practitioners must evaluate whether the timescale of the target biophysical transition (e.g., sidechain reorientations requiring microseconds, or domain opening taking milliseconds) is accessible with available compute before committing to atomistic MD (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

## Thermodynamic ensembles, thermostats, and barostats

Simulations must strictly generate the correct statistical mechanical probability distribution for the chosen ensemble:
- **Microcanonical ($NVE$)**: Represents an isolated system with constant particle number, volume, and total energy; used primarily for testing integrator stability and verifying energy conservation (absence of drift) (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Canonical ($NVT$)**: Constant particle number, volume, and temperature; essential for gas-phase, surface, or fixed-density systems, and as the initial thermalization phase for condensed systems (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Isothermal–Isobaric ($NPT$)**: Constant particle number, pressure, and temperature; the standard production ensemble for condensed biomolecular systems matching laboratory bench conditions (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

### Thermostat Audit and Known Artifacts

Braun et al. categorize and evaluate popular temperature-regulation algorithms:

| Thermostat | Ensemble | Nature | Space | Physical? | Structural Properties | Dynamical Properties | Known Shortcomings / Artifacts |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **None** | $NVE$ | Deterministic | Global | Yes | Correct | Correct | Unsuitable for temperature regulation |
| **Simple Velocity Rescaling** | Undefined | Deterministic | Global | No | Incorrect | Incorrect | Suppresses kinetic energy fluctuations |
| **Berendsen** | Undefined | Deterministic | Global | No | Incorrect | Incorrect | Violates equipartition theorem; prone to the **"flying ice cube"** artifact where kinetic energy drains from high-frequency modes into center-of-mass translation |
| **Bussi (Stochastic Velocity Rescaling / v-rescale)** | Canonical ($NVT$) | Stochastic | Global | No | Correct | Mostly correct | Introduces a stochastic Wiener process into velocity scaling, ensuring the exact canonical ensemble while preserving momentum and physical dynamics |
| **Andersen** | Canonical ($NVT$) | Stochastic | Local | No | Correct | Incorrect | Random velocity re-assignments disrupt velocities and momentum autocorrelation functions; distorts transport properties (diffusion coefficients) |
| **Langevin Dynamics** | Canonical ($NVT$) | Stochastic | Local | No | Correct | Incorrect | Introduces friction and random noise terms; alters physical time-correlation functions unless friction coefficient $\gamma$ is very small |
| **Nosé–Hoover** | Canonical ($NVT$) | Deterministic | Global | Yes | Correct | Mostly correct | Extended Hamiltonian with a fictitious thermal reservoir; can exhibit non-ergodic sampling for stiff or small systems unless Nosé–Hoover chains (NHC) are implemented |

(source: Best Practices for Foundations in Molecular Simulations [Article.pdf)

### Barostat Selection and Pressure Regulation

- **Berendsen Barostat**: Scales coordinates with an empirical relaxation time $\tau_p$; does not sample the true isothermal–isobaric ensemble and suppresses volume fluctuations (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). It should be restricted to early equilibration stages and strictly avoided for production simulations or property calculations sensitive to compressibility (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Parrinello–Rahman and Andersen Barostats**: Extended-ensemble formulations that yield the exact $NPT$ ensemble; prone to numerical instabilities or box collapse if launched from poorly equilibrated systems with high initial stress (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Monte Carlo Barostat**: Attempts periodic trial volume changes accepted via the Metropolis criterion; efficiently yields the true $NPT$ ensemble without introducing artificial dynamical coupling into the equations of motion (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

## Integration timesteps and boundary conditions

- **Timestep Selection Rule**: The integration timestep $\Delta t$ should not exceed approximately one-tenth of the period of the fastest motion in the system ($\Delta t \le \tau_{\text{fastest}} / 10$) (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). In all-atom systems, bond stretching involving hydrogen atoms occurs on a $\sim 10\text{ fs}$ timescale, setting a limit of $\Delta t \le 1.0\text{ fs}$ for fully flexible models (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). Applying holonomic constraints (SHAKE, RATTLE, or LINCS) to all hydrogen bonds removes these fast motions and allows $\Delta t = 2.0\text{ fs}$ (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Energy Drift Diagnostics**: The chosen timestep must be validated by running a short simulation in the microcanonical ($NVE$) ensemble; any noticeable upward or downward drift in total energy denotes integrator breakdown or inadequate cutoff settings (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Periodic Boundary Conditions (PBC)**: Simulated cells represent an infinite lattice under minimum image convention (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). To prevent a solute from interacting with its own periodic image, the shortest box dimension must be at least twice the nonbonded cutoff distance ($L_{\text{box}} \ge 2 r_{\text{cut}}$) (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Cutoffs and Long-Range Forces**: Electrostatic interactions are long-range ($r^{-1}$) and must be treated with lattice-sum techniques (Particle Mesh Ewald, PME) rather than simple truncation (source: Best Practices for Foundations in Molecular Simulations [Article.pdf). Van der Waals interactions ($r^{-6}$) typically employ cutoffs between $1.0\text{ and }1.5\text{ nm}$, and analytical dispersion corrections should be applied to account for neglected long-range attraction when computing densities and surface tensions (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

## Equilibration, convergence, and replicate sampling

Braun et al. emphasize that thermodynamic equilibration cannot be inferred from a single observable:
- Energy, temperature, and pressure may equilibrate rapidly (within picoseconds to nanoseconds), whereas structural relaxation (RMSD, pocket opening, sidechain rotamers, lipid bilayer packing) can take hundreds of nanoseconds (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- Discarding initial non-equilibrium data ("burn-in" period) is mandatory (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
- **Independent Replicates**: Running multiple independent simulations initiated with different random velocity seeds from the same or distinct starting configurations is far more effective at diagnosing sampling convergence and estimating robust statistical uncertainties than relying on a single extended trajectory (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

## The LiveCoMS foundational checklist

The paper synthesizes a checklist for simulation planning:
1. **Take stock of plans**: Define reference states, target experimental observables, temperature, pressure, force field suitability, and timescale feasibility (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
2. **System preparation**: Verify chemical components, protonation states, force field assignments, and box dimensions (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
3. **Cutoff handling**: Ensure $L_{\text{box}} \ge 2 r_{\text{cut}}$, choose appropriate PME parameters, and enable long-range dispersion corrections (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
4. **Ensemble settings**: Choose thermostats and barostats that reproduce the true thermodynamic ensemble (avoiding Berendsen in production runs) (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
5. **Timestep stability**: Verify timestep choice via $NVE$ energy conservation checks (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).
6. **Run protocol & storage**: Establish phased minimization ($NVT \to NPT \to \text{production}$), determine frame saving frequency to balance statistical independence with storage capacity, and execute replicate trajectories (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

## Related pages
- [[simulation-and-benchmark-best-practices]]
- [[molecular-dynamics-simulation]]
- [[force-fields]]
- [[gromacs-protein-ligand-protocol]]
- [[lemkul-gromacs-protein-ligand-tutorial]]
- [[hydrogen-mass-repartitioning]]
- [[binding-free-energy]]
- [[hahn-2022-binding-affinity-benchmarks]]

