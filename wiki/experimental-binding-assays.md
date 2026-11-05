# Experimental Binding Assays

**Summary**: A guide to biophysical binding assays—chiefly Isothermal Titration Calorimetry and Surface Plasmon Resonance—used to measure thermodynamic and kinetic binding parameters and validate computational predictions.
**Sources**: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md; The Light and Dark Sides of Virtual Screening.pdf; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf
**Last updated**: 2026-11-01

---

## Role in computational drug discovery

Experimental biophysical binding assays provide the ground-truth benchmarks required to calibrate [[scoring-functions]], validate [[protein-ligand-docking]] poses, and assess the accuracy of [[binding-free-energy]] and [[binding-kinetics]] calculations (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; The Light and Dark Sides of Virtual Screening.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md). While in silico methods estimate binding thermodynamics and kinetic rates, experimental measurements confirm physical target engagement, determine stoichiometry, calibrate standard-state volume corrections ($V_0 = 1661\text{ \AA}^3$), and eliminate false-positive computational artifacts (source: The Light and Dark Sides of Virtual Screening.pdf; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

The two primary quantitative techniques in structure-based drug design are **Isothermal Titration Calorimetry (ITC)** and **Surface Plasmon Resonance (SPR)** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Thermodynamics and Kinetics of Drug-Target Binding by Molecular Simulation.md).

## Isothermal titration calorimetry (ITC)

ITC is the gold-standard technique for measuring binding thermodynamics in solution without requiring chemical labeling or surface immobilization (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### 1. Instrument design and operational principle
An ITC instrument contains a reference cell filled with buffer or water and an adiabatic sample cell containing the macromolecular receptor at concentration $[M]_t$ (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). A computer-controlled syringe delivers a series of precise aliquots (typically 1–10 $\mu$L) of the ligand solution into the stirred sample cell:
- **Heat of reaction**: Each injection generates or absorbs a burst of heat ($q_i$) proportional to the amount of ligand bound (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Thermal feedback**: A feedback power heater maintains the sample cell at the exact temperature of the reference cell. The instrumental signal is the electric power (in $\mu\text{cal/s}$ or $\mu\text{W}$) required to null the temperature difference between cells (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Thermogram integration**: Integrating the power peaks over time yields the heat per injection, which decreases as the receptor approaches saturation (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### 2. Extracted thermodynamic parameters
A single ITC titration experiment directly resolves three parameters by non-linear regression fitting to the binding isotherm (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):
1. **Binding constant ($K_b$ or $K_a = 1/K_d$)**;
2. **Enthalpy of binding ($\Delta H^\circ$)**;
3. **Binding stoichiometry ($n$)** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

From these values, the Gibbs free energy and conformational/solvation entropy are directly computed:
$$\Delta G^\circ = -RT \ln K_b$$
$$-T\Delta S^\circ = \Delta G^\circ - \Delta H^\circ$$
(source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). Furthermore, executing titrations across multiple temperatures yields the binding heat capacity change:
$$\Delta C_p = \frac{\partial \Delta H}{\partial T}$$
Negative $\Delta C_p$ values provide a signature of the hydrophobic effect and non-polar surface burial (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). See [[binding-thermodynamics]] and [[water-and-solvation]].

### 3. The Wiseman parameter ($c$)
Accurate determination of $K_b$ requires that the dimensionless Wiseman parameter $c$ falls within an optimal experimental window (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md):
$$c = n K_b [M]_t = \frac{n [M]_t}{K_d}$$
- For reliable deconvolution of $K_b$ and $\Delta H$, $c$ should lie between **$10$ and $500$** (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- If $c > 1,000$ (very tight binding, $K_d < 1\text{ nM}$), the titration isotherm becomes a step function, accurately measuring $\Delta H$ and $n$ but preventing extraction of $K_b$. Displacement titrations using a competing weak ligand are required in this regime (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- If $c < 1$ (weak binding, $K_d > 100\text{ }\mu\text{M}$), the curve lacks an inflection point, complicating non-linear fitting (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### 4. Technical limitations of ITC
- **High sample consumption**: Demands substantial quantities of purified protein (often milligrams per run at concentrations of 10–100 $\mu$M) (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Throughput constraints**: Runs are serial and require thorough cell cleaning, limiting throughput to 2–4 titrations per day on standard instruments (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Enthalpy of dilution & DMSO sensitivity**: Heats of ligand dilution and DMSO mismatch between syringe and cell must be subtracted to prevent baseline distortions (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

## Surface plasmon resonance (SPR)

SPR is an optical biosensing technology that monitors real-time macromolecular interactions to extract both kinetic rate constants ($k_{\text{on}}, k_{\text{off}}$) and thermodynamic constants ($K_d$) (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).

### 1. Physical foundation and sensorgram
- **Kretschmann optical configuration**: Polarized light illuminates the back of a sensor chip coated with a thin gold film (typically ~50 nm). Under total internal reflection, an evanescent wave excites collective surface electron oscillations (plasmons) in the gold layer, producing a sharp dip in reflected light intensity at the resonance angle $\theta_{\text{SPR}}$ (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Sensor response**: The resonance angle depends sensitively on the local refractive index within ~200–300 nm of the sensor surface. When target protein is immobilized on the sensor dextran matrix and ligand flows through the microfluidic channel, binding increases local mass and refractive index, shifting $\theta_{\text{SPR}}$ proportionally in Resonance Units (RU, where $1\text{ RU} \approx 1\text{ pg/mm}^2$) (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Sensorgram phases**:
  1. *Baseline*: Continuous buffer flow establishes the initial signal;
  2. *Association phase*: Ligand solution flows over the sensor surface; binding produces a rising curve governed by $k_{\text{on}}[L] + k_{\text{off}}$;
  3. *Steady state*: Equilibrated plateau where association rate equals dissociation rate;
  4. *Dissociation phase*: Buffer replaces the ligand solution, and the signal decays exponentially according to $k_{\text{off}}$, enabling direct calculation of residence time $\tau = 1/k_{\text{off}}$ (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf).

### 2. Kinetic and affinity deconvolution
By fitting sensorgrams across a concentration series to a 1:1 Langmuir binding model, SPR resolves:
- Association rate constant ($k_{\text{on}}$, in $\text{M}^{-1}\text{s}^{-1}$);
- Dissociation rate constant ($k_{\text{off}}$, in $\text{s}^{-1}$);
- Equilibrium dissociation constant:
  $$K_d = \frac{k_{\text{off}}}{k_{\text{on}}}$$
  (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md; Understanding the impact of binding free energy and kinetics calculations in modern drug discovery.pdf). See [[binding-kinetics]].

### 3. Experimental artifacts and pitfalls in SPR
- **Mass transport limitation (MTL)**: When the intrinsic association rate is faster than the diffusion of ligand from bulk flow to the sensor surface ($k_{\text{on}} > 10^6\text{ M}^{-1}\text{s}^{-1}$), the observed rate reflects diffusion rather than chemistry, distorting kinetic fitting unless flow rates are elevated and immobilization levels kept low (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Surface immobilization artifacts**: Covalent amine coupling (EDC/NHS) or sensor tethering can alter protein conformational dynamics, block binding epitopes, or restrict rotational entropy, altering measured rates relative to solution experiments (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Bulk refractive index shifts**: Differences in DMSO concentration (e.g., $\Delta \text{DMSO} \approx 0.1\%$) between running buffer and sample produce large step shifts in RU that require solvent calibration curves (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md).
- **Non-specific binding and colloidal aggregation**: Promiscuous aggregators form sub-micrometer colloidal particles that adhere non-specifically to sensor surfaces. Including non-ionic detergents (e.g., 0.01%–0.05% P20 or Triton X-100) and running control reference surfaces are required to identify false positives (source: The Light and Dark Sides of Virtual Screening.pdf). See [[docking-failure-modes]].

## Comparison of biophysical methods

| Technique | Primary Readouts | Protein Amount | Throughput | Immobilization? | Strengths & Vulnerabilities |
|---|---|---|---|---|---|
| **ITC** | $\Delta H, K_b, n, \Delta G, -T\Delta S, \Delta C_p$ | High (mg) | Low (serial) | No (solution) | Direct thermodynamic decomposition; high protein consumption (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). |
| **SPR** | $k_{\text{on}}, k_{\text{off}}, K_d, \tau$ | Low ($\mu$g) | Medium-High | Yes (chip) | Direct kinetic rate profiling; subject to mass transport and immobilization artifacts (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). |
| **NMR (CSP)** | Residue-level $K_d$, pocket localization | Medium-High | Medium | No (solution) | Identifies binding site at atomic resolution; limited by protein molecular weight (source: Insights into Protein–Ligand Interactions_ Mechanisms, Models, and Methods.md). |
| **DSF (Thermal Shift)** | $\Delta T_m$ (melting shift) | Very Low ($\mu$g) | Very High (384-well) | No (solution) | Fast screening for hit identification; semi-quantitative affinity estimation (source: The Light and Dark Sides of Virtual Screening.pdf). |

## Related pages
- [[binding-thermodynamics]]
- [[binding-kinetics]]
- [[binding-free-energy]]
- [[structure-based-drug-design]]
- [[virtual-screening]]
- [[docking-failure-modes]]
- [[protein-ligand-interactions]]
- [[chemical-screening-artifacts]]
