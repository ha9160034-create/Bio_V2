# Enhanced Sampling

**Summary**: Enhanced sampling methods accelerate molecular dynamics so that rare events — ligand binding and unbinding, pocket opening, conformational transitions — occur within accessible simulation time, usually by adding biases, exchanging replicas, scaling the potential, or running many short trajectories.
**Sources**: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulations and drug discovery.md; Molecular dynamics simulation for all.md; Enhanced Sampling Techniques in Molecular Dynamics Simulations of Biological Systems.pdf; Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf; Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf
**Last updated**: 2026-11-01

---

## The sampling problem

MD with a 1 fs time step requires 10¹² integration steps to cover a millisecond and 10¹⁵ to cover a second. Ligand binding and unbinding typically require milliseconds to seconds, far beyond routine simulation. GPUs make microseconds routine, and Anton reaches milliseconds, but the second timescale remains elusive (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; Molecular dynamics simulations and drug discovery.md). Furthermore, static crystal structures capture only deep local energy minima; discovering cryptic pockets or induced-fit conformations requires crossing substantial free-energy barriers that standard microsecond MD rarely traverses (source: Molecular dynamics simulations and drug discovery.md).

Two distinct needs must be met (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md):

1. **Sample a complex landscape** of hierarchically organized basins separated by barriers, where ergodicity is nearly broken.
2. **Quantify kinetics** under stationary or nonequilibrium conditions — harder than equilibrium sampling.

Enhanced sampling methods share a way of accelerating events of interest, often by focusing sampling in specific regions using **collective variables** (CVs): distances, angles, RMSD, or other observables ξ whose value tracks the process (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Steered MD

A parabolic potential is added to the plain MD potential, and its center moves along the reaction coordinate at constant velocity:

$$U_{bias} = \frac{k}{2}[\xi - \xi_0(t)]^2, \qquad \xi_0(t) = vt$$

The ligand is effectively pulled out of the site. Park and Schulten showed how to recover the potential of mean force (PMF) from such nonequilibrium simulations using Jarzynski's equality, by running several independent replicas of the same steering process (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Steered MD gives a free-energy profile along a chosen path and can reveal barriers, but the result depends on the pulling speed and the choice of coordinate.

## Adiabatic bias MD (ABMD)

ABMD applies a harmonic restraint whose center moves only when the collective variable advances in the desired direction — a pawl-and-ratchet mechanism. If the system tries to move backward, the restraint prevents it. Unlike steered MD, the speed is set partly by the system's natural evolution; the restraint constant tunes how fast the center advances. ABMD can be viewed as an adaptive, gentle version of steered MD (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Parallel tempering

Identical replicas are run at different temperatures, with the base temperature T₁ and higher temperatures T_i. Configurations are exchanged between adjacent replicas with a Metropolis acceptance probability, allowing high-temperature configurations to migrate into the low-temperature replica and escape basins that the low-temperature simulation alone would not visit (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Equal acceptance ratios require a geometric progression of temperatures. Parallel tempering satisfies detailed balance. Replica exchange can also swap other observables, such as metadynamics bias potentials (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Scaled MD

A positive constant μ ∈ (0, 1] defines a new potential:

$$U_\mu(\mathbf{r}) = \mu U(\mathbf{r})$$

Scaling uniformly flattens all barriers; μ = 1 is plain MD, and smaller μ gives more smoothing. A constant scaling is equivalent to a 1/μ scaling of temperature (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

For unbinding, harmonic restraints are applied to the protein backbone outside the binding site to prevent unfolding. The time at which the ligand becomes fully surrounded by water is defined as the unbinding time; at least 20 repeats are used to define an average. Compounds are ranked by this time. The method needs few free parameters and no reaction coordinate, but it requires prior knowledge of the binding-site residues, and heavy scaling can distort the unbinding pathway (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

A published k_off ranking used λ = 0.4 on HSP90, Grp78, and the A2A GPCR with three ligand series each (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). See [[binding-kinetics]].

## τ-RAMD

Random acceleration MD (RAMD) periodically applies a random force to the ligand. If the ligand does not move in the desired direction — judged by a distance threshold — the force direction is reassigned. This accelerates unbinding by orders of magnitude. τ-RAMD runs, on average, **40–200 simulations per compound** and uses the mean unbinding time τ to build correlations with experimental residence times (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

It is a supervised method in the sense that the random force is coupled to an obvious progress coordinate, the ligand–site distance.

## Metadynamics

Metadynamics adds a history-dependent bias built from Gaussians deposited along the collective variables:

$$V(\xi,t) = \sum_{t'=\tau_G, 2\tau_G,\dots < t} W \exp\left[-\sum_{i=1}^{d}\frac{(\xi_i-\xi_i(t'))^2}{2\sigma_i^2}\right]$$

where ω = W/τ_G is the energy rate, σ_i is the Gaussian width, and τ_G is the deposition stride. The bias fills basins sequentially; at convergence, the negative of the bias approximates the free-energy surface, F(ξ) ≈ −V(ξ) + C (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Advantages: it explores CV space automatically and yields the free energy at convergence. Drawbacks: it is hard to know when to stop, the last deposited Gaussian introduces an arbitrary offset, and there is no fully satisfying convergence proof. **Well-tempered metadynamics** fixes the stopping problem by reducing the Gaussian height over time and has been proved to converge (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Metadynamics has been applied to ligand unbinding, residence-time estimation, and binding-site discovery; see [[binding-kinetics]].

## Umbrella sampling

Umbrella sampling samples with a harmonic restraint centered at successive positions along the reaction coordinate:

$$U_{biased} = U(\mathbf{r}) + \frac{k}{2}[\xi(\mathbf{r})-\xi_0]^2$$

Many windows are run to cover the coordinate range, and the windows are recombined into a single PMF. **WHAM** (Weighted Histogram Analysis Method) and **MBAR** (Multistate Bennett Acceptance Ratio) are the standard recombinations. A second class of methods extracts the mean force directly, adapting thermodynamic integration to umbrella sampling (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

To compare with experiment, the PMF must be converted to a standard binding free energy. The standard-state volume V₀ = 1661 Å³ corresponds to C₀ = 1 M; a correction is required to move from that volume to the unbound volume actually sampled. Restraints orthogonal to the reaction coordinate add another free-energy term (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Umbrella sampling is well controlled but requires many windows and a good reaction coordinate.

## Adaptive biasing force (ABF)

ABF estimates the mean force on the reaction coordinate on the fly and applies a bias opposing it, so the system can escape local minima. As the running average converges to the true mean force, the total force felt by the system vanishes and the coordinate diffuses over its range; the free energy is then obtained by integrating the mean force (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Markov state models (MSM)

An MSM partitions configuration space into non-overlapping states and estimates a transition probability matrix T_ij over a lag time τ. The lag time must be long enough for transitions to be memoryless and short enough for resolution; the Chapman–Kolmogorov test validates the choice (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Rows sum to one. The model predicts long-term evolution from an out-of-equilibrium state and yields equilibrium populations, free energies, and rates. MSMs work best for diffusive problems and less well for high-barrier activated events, because plain MD may never cross the barrier. Refinement by targeted short runs is often needed (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Examples include 495 × 100 ns trypsin–benzamidine simulations (37% spontaneous binding) and a 150 µs MSM that identified six apo-trypsin states with different surface pockets (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Weighted ensemble (WE) and WExplore

The weighted ensemble method partitions the path from state A to state B into bins. It starts M trajectories in the bin containing A, advances them for a fixed short time τ, then prunes and spawns trajectories to keep the number of walkers per bin constant. The cycle repeats until B is reached. The underlying dynamics are unbiased (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

WE provides unbiased predictions on timescales longer than the aggregated simulation time, based on the Hill relation between mean first-passage time and steady-state flux. It handles multiple pathways, parallelizes easily, and does not require a Markovian assumption (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

**WExplore** defines bins dynamically and hierarchically as Voronoi polyhedra, avoiding an a priori bin definition and reducing a high-dimensional order-parameter space to a manageable tree (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

Notable results:

- TPPU unbinding from soluble epoxide hydrolase: experimental residence time 11 min (660 s); WExplore estimated 42 s with a standard error of ~10² s from only 6 µs of aggregated simulation — within 1–2 orders of magnitude (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).
- Benzamidine from trypsin: predicted exit time 180 µs versus measured 1700 µs, and revealed three exit channels, one previously unknown (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Classification: CV-based versus non-CV-based methods

Enhanced sampling methodologies are broadly organized into two design paradigms (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf):
- **Collective variable (CV)-based**: Methods such as Umbrella Sampling, Adaptive Biasing Force (ABF), Steered MD, and Metadynamics require user-defined reaction coordinates (distances, angles, RMSD).
- **Non-CV-based (reaction-coordinate-free)**: Methods such as Replica Exchange (REMD), Accelerated MD (aMD), and Gaussian Accelerated MD (GaMD) accelerate sampling across the entire potential energy surface without requiring prior knowledge of transition pathways (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).

### Gaussian accelerated MD (GaMD) and LiGaMD
Accelerated MD (aMD) flattens energy barriers by adding a non-negative boost potential $V_{\text{boost}}$ when potential energy $V < E$, but large potential fluctuations hinder accurate free energy reweighting (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).
- **GaMD formulation**: Calibrates harmonic boost potentials so that the boosted potential energy follows a near-Gaussian distribution (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).
- **Cumulant expansion reweighting**: The Gaussian distribution enables unbiasing using cumulant expansion to the second order, recovering unperturbed free energy surfaces without reweighting noise (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf).
- **Ligand GaMD (LiGaMD)**: Applies dual boost potentials to the ligand–protein non-bonded interaction and ligand bonded dihedrals, capturing repetitive spontaneous binding and unbinding on microsecond timescales to calculate both $\Delta G_{\text{bind}}$ and kinetic rates ($k_{\text{on}}, k_{\text{off}}$) (source: Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf; Enhanced-Sampling Simulations for the Estimation of Ligand Binding Kinetics Current Status and Perspective.pdf).

### Hydrogen mass repartitioning (HMR)
HMR transfers mass from bonded heavy atoms (C, N, O) to attached hydrogen atoms, increasing hydrogen mass to ~3.0 amu (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Molecular dynamics simulation for all.md):
- Slows high-frequency C–H and O–H stretching vibrations;
- Permits expanding the numerical integration timestep from 2.0 fs to **4.0 fs** in Verlet integrators when combined with holonomic bond constraints (SHAKE/LINCS), yielding a **2-fold speedup** in simulation throughput without altering conformational free energy landscapes or thermodynamic order parameters (source: Advances in Molecular Dynamics Simulations and Enhanced Sampling Methods for the Study of Protein Systems.pdf; Molecular dynamics simulation for all.md). See [[hydrogen-mass-repartitioning]].

### Generalized simulated annealing (GSA)
Developed within Tsallis non-extensive thermostatistics, GSA is governed by visitation ($q_v$) and acceptance ($q_A$) indices (source: Enhanced Sampling Techniques in Molecular Dynamics Simulations of Biological Systems.pdf):
- Facilitates frequent local sampling coupled with occasional long-distance exploratory jumps in high-dimensional phase space;
- Outperforms classical simulated annealing (CSA) and fast simulated annealing (FSA) in escaping deep energy minima during peptide folding and force-field parameterization (source: Enhanced Sampling Techniques in Molecular Dynamics Simulations of Biological Systems.pdf).

### Supervised molecular dynamics (SuMD)
Monitors ligand–cavity distances across successive short classical MD windows (e.g., 600 ps) (source: Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf):
- Advances the simulation only when the distance decreases, reassigning velocities if progress stalls;
- Tracks complete unbinding-to-binding pathways and intermediate desolvation steps in nanoseconds without altering the potential energy function (source: Bridging Molecular Docking to Molecular Dynamics in Exploring Ligand-Protein Recognition Process An Overview.pdf).

## Choosing a method

| Goal | Reasonable choice |
|---|---|
| Free-energy profile along a known coordinate | Umbrella sampling, ABF, metadynamics |
| Pulling/unbinding pathway | Steered MD, ABMD, τ-RAMD |
| Residence-time ranking without a coordinate | Scaled MD, τ-RAMD, WExplore |
| Equilibrium populations and rates from many short runs | MSM |
| Long-timescale kinetics with unbiased dynamics | Weighted ensemble / WExplore |
| Enhanced configurational sampling without a CV | Parallel tempering, accelerated MD |

Guidance from the literature: methods whose outcome depends heavily on many human choices are not mature enough for prospective quantitative use; independent runs are needed to assess reliability; and the choice is system-dependent (source: Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Related pages
- [[molecular-dynamics-simulation]]
- [[hydrogen-mass-repartitioning]]
- [[binding-free-energy]]
- [[alchemical-free-energy-methods]]
- [[mm-pbsa-and-mm-gbsa]]
- [[binding-kinetics]]
- [[force-fields]]
- [[cryptic-and-allosteric-pockets]]
- [[benchmark-drug-targets]]
- [[decherchi-2020-binding-thermodynamics-kinetics]]
- [[bernardi-2015-enhanced-sampling-md]]
- [[lazim-2020-md-enhanced-sampling-advances]]
- [[ahmad-2022-enhanced-sampling-binding-kinetics]]
- [[metadynamics]]
- [[umbrella-sampling]]
- [[markov-state-models]]
