# GROMACS Protein–Ligand Simulation Protocol

**Summary**: A comprehensive, step-by-step technical protocol for setting up, equilibrating, running, and analyzing all-atom protein–ligand molecular dynamics simulations in GROMACS, preserving exact command syntax, parameterization rules, and analysis workflows.
**Sources**: `gromacs_protein_ligand_complex_tutorial_lemkul.md`, `mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf`, `Best Practices for Foundations in Molecular Simulations [Article.pdf`
**Last updated**: 2026-11-02

---

## Methodological framework

Simulating a protein–ligand complex in classical molecular dynamics requires decoupling the macromolecular receptor from the small-molecule ligand during initial topology generation (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). Because standard force-field databases (`.rtp`) contain parameter definitions strictly for canonical amino acids and nucleic acids, small molecules require external parameterization prior to system recombination (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). 

This protocol synthesizes the standardized workflow formulated by Justin Lemkul and the LiveCoMS best-practice guidelines for GROMACS, utilizing the **CHARMM36 all-atom force field** coupled to the **CHARMM General Force Field (CGenFF)** as the primary reference, while providing cross-compatibility notes for AMBER/GAFF workflows (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).

## Step 1: Receptor preparation and protein topology (`pdb2gmx`)

1. Obtain the experimental crystallographic complex (e.g., T4 lysozyme L99A/M102Q bound to 2-propylphenol, PDB ID `3HTB`) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
2. Clean the coordinate file using a molecular viewer (PyMOL, Chimera) or text parsing: strip non-structural crystallization additives (e.g., phosphate `PO4`, $\beta$-mercaptoethanol `BME`) and non-bridging bulk crystallographic waters (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). *(Note: Structurally conserved water molecules bridging ligand–protein hydrogen bonds must be preserved; see [[receptor-preparation]])* (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
3. Extract the ligand coordinates into a separate PDB file and remove the ligand records from the protein file:
   ```bash
   grep "JZ4" 3HTB_clean.pdb > jz4.pdb
   grep -v "JZ4" 3HTB_clean.pdb > 3HTB_protein.pdb
   ```
   (source: gromacs_protein_ligand_complex_tutorial_lemkul.md)
4. Unarchive the latest CHARMM36 force field directory in the working directory (e.g., `charmm36-feb2026_cgenff-5.0.ff`) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
5. Generate the protein topology using `gmx pdb2gmx`:
   ```bash
   gmx pdb2gmx -f 3HTB_protein.pdb -o 3HTB_processed.gro -p topol.top -i posre.itp -ter
   ```
   - Select **CHARMM36 all-atom force field** (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
   - Select the CHARMM-modified **TIP3P** water model (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
   - Assign terminal capping patches (`NH3+` and `COO-` for zwitterionic termini) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
   - Note the reported total charge of the macromolecule (e.g., net charge **+6e**) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Step 2: Small-molecule ligand parameterization

### 1. Preparation of `.mol2` Input
CHARMM is an all-atom force field requiring explicit hydrogen atoms assigned to correct valences and hybridization (source: gromacs_protein_ligand_complex_tutorial_lemkul.md):
- Open `jz4.pdb` in Avogadro, execute **Build $\to$ Add Hydrogens**, and save as Sybyl Mol2 (`jz4.mol2`) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Edit `jz4.mol2` in a text editor: under `@<TRIPOS>MOLECULE`, replace `*****` with the residue identifier `JZ4`; in `@<TRIPOS>ATOM`, unify all atom lines so that both heavy atoms and hydrogens belong to residue `1 JZ4` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- Sort the bond section into ascending order using the helper script:
  ```bash
  perl sort_mol2_bonds.pl jz4.mol2 jz4_fix.mol2
  ```
  (source: gromacs_protein_ligand_complex_tutorial_lemkul.md)

### 2. CGenFF Submission and Penalty Score Validation Rules
Submit `jz4_fix.mol2` to the CGenFF server (`cgenff.com` or `paramchem.org`) to generate the CHARMM stream file (`jz4.str`) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). Inspect the penalty scores assigned to charges and dihedral angles:
- **Penalty $< 10$**: Parameters are robust, transferable, and ready for production MD without modification (for JZ4, penalties are $< 1$) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- **Penalty $10\text{ to }50$**: Parameters require manual verification, such as QM water interaction energy calculations or dihedral potential energy scans (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
- **Penalty $> 50$**: Parameters are highly dubious; extensive quantum chemical (DFT) re-parameterization of partial charges and torsional profiles is mandatory (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

*(Alternative AMBER/GAFF Workflow: When simulating under AMBER force fields such as Amber14SB, use Antechamber with AM1-BCC partial charges and `acpype` to convert GAFF/GAFF2 parameters into GROMACS `.itp` format: `acpype -i jz4.mol2 -c bcc -b JZ4`)* (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

### 3. Converting and Structuring the Ligand `.itp` File
Convert the CGenFF stream file into GROMACS format (`jz4_fix_gmx.top` and `jz4_fix_gmx.pdb`) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). Copy `jz4_fix_gmx.top` to `jz4_gmx.itp` and format it as an included topology module:
1. Strip `#include "charmm36.ff/forcefield.itp"` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
2. Paste the `[ dihedraltypes ]` contents directly above `[ moleculetype ]` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
3. Set the `[ moleculetype ]` name to `JZ4` with `nrexcl = 3` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
4. Point position restraints to `#include "posre_jz4.itp"` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
5. Strip all solvent and system directives from the bottom of the file (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

### 4. Merging Coordinates and Master Topology
1. Convert ligand coordinates to GROMACS `.gro` format:
   ```bash
   gmx editconf -f jz4_fix_gmx.pdb -o jz4.gro
   ```
   (source: gromacs_protein_ligand_complex_tutorial_lemkul.md)
2. Copy `3HTB_processed.gro` to `complex.gro`. Append the 22 atom lines from `jz4.gro` immediately below the final protein residue and update the total atom count on line 2 ($2614 + 22 = 2636$ atoms) (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
3. In `topol.top`, insert `#include "jz4_gmx.itp"` immediately after the parent force field `#include` and before the protein `[ moleculetype ]` directive (source: gromacs_protein_ligand_complex_tutorial_lemkul.md):
   ```ini
   ; Include forcefield parameters
   #include "./charmm36-feb2026_cgenff-5.0.ff/forcefield.itp"

   ; Include ligand parameters and topology
   #include "jz4_gmx.itp"

   [ moleculetype ]
   ; Name            nrexcl
   Protein_chain_A     3
   ```
   (source: gromacs_protein_ligand_complex_tutorial_lemkul.md)
4. Append `JZ4 1` to the bottom of `topol.top` under `[ molecules ]` (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Step 3: Box definition, solvation, and charge neutralization

1. **Unit Cell Geometry**: Define a **rhombic dodecahedral box** with a $1.0\text{ nm}$ buffer distance between solute and box boundaries (`-bt dodecahedron -d 1.0`):
   ```bash
   gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
   ```
   *Physical Justification*: Rhombic dodecahedral geometry occupies only $\sim 77\%$ of the volume of a cubic box with identical minimum periodic image distance, reducing explicit solvent waters by $\sim 23\%$ and accelerating calculation speed proportionally (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; source: mshirts,+5068-from-proteins-to-perturbed-hamiltonians-a-suite-of-tutorials-for-the-gromacs-2018-molecular-simulation-package-article.pdf).
2. **Solvation**: Fill the unit cell with explicit water molecules:
   ```bash
   gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
   ```
   (source: gromacs_protein_ligand_complex_tutorial_lemkul.md)
3. **Ion Neutralization**: Generate the run input for ion placement and replace solvent molecules with counterions to reach net electroneutrality (mandatory for Particle Mesh Ewald electrostatics):
   ```bash
   gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
   gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
   ```
   Select group `SOL` when prompted (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Step 4: Energy minimization

Eliminate initial steric clashes and severe geometric strains using steepest descents (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

### Force-Switch Nonbonded Settings for CHARMM36 (`em.mdp`)
```ini
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
nstlist         = 20
rlist           = 1.2
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw            = 1.2
coulombtype     = PME
rcoulomb        = 1.2
pme_order       = 4
fourierspacing  = 0.16
```
*(AMBER Setting Note: Under AMBER force fields, replace `vdw-modifier = force-switch` with standard cutoff `vdw-modifier = none`, set `rvdw = 0.9` or `1.0` nm, and enable analytical dispersion corrections `DispCorr = EnerPres`)* (source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

Execute minimization:
```bash
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```
Verify convergence: maximum force $F_{\text{max}} < 1000\text{ kJ/(mol}\cdot\text{nm)}$ (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

## Step 5: Equilibration with position restraints

### 1. Generating Restraints for Ligand Heavy Atoms
To prevent the unbonded ligand from drifting out of the binding pocket while solvent and ions relax, apply harmonic position restraints ($k = 1000\text{ kJ/(mol}\cdot\text{nm}^2)$) to ligand heavy atoms (source: gromacs_protein_ligand_complex_tutorial_lemkul.md):
```bash
gmx make_ndx -f jz4.gro -o index_jz4.ndx
# In prompt: 0 & ! a H* -> name group "JZ4_Heavy" -> q
gmx genrestr -f jz4.gro -n index_jz4.ndx -o posre_jz4.itp -fc 1000 1000 1000
```
(source: gromacs_protein_ligand_complex_tutorial_lemkul.md)

### 2. Merging Temperature Coupling Groups (`tc-grps`)
A critical best-practice rule: **never couple a small ligand to its own thermostat or lump it with bulk solvent** (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; source: Best Practices for Foundations in Molecular Simulations [Article.pdf). The 22 atoms of JZ4 lack sufficient degrees of freedom for independent thermostatting, causing severe kinetic temperature fluctuations; conversely, grouping it with `non-Protein` distorts solvent coupling (source: gromacs_protein_ligand_complex_tutorial_lemkul.md). Create a merged `Protein_JZ4` index group:
```bash
gmx make_ndx -f em.gro -o index.ndx
# In prompt: 1 | 13 (Protein | JZ4) -> name group "Protein_JZ4" -> q
```
In both `nvt.mdp` and `npt.mdp`, configure:
```ini
tc-grps = Protein_JZ4 Water_and_ions
tau_t   = 0.1         0.1
ref_t   = 300         300
```
(source: gromacs_protein_ligand_complex_tutorial_lemkul.md)

### 3. NVT and NPT Equilibration Runs
Execute 100 ps of $NVT$ thermalization (stabilizing temperature at $300\text{ K}$ using the Bussi stochastic velocity-rescaling thermostat `tcoupl = v-rescale`), followed by 100 ps of $NPT$ isobaric equilibration (stabilizing pressure at $1.0\text{ bar}$ and solvent density at $\sim 1000\text{ kg/m}^3$ using the Berendsen barostat `pcoupl = berendsen`):
```bash
# NVT Equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT Equilibration
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
```
(source: gromacs_protein_ligand_complex_tutorial_lemkul.md)

## Step 6: Unrestrained production molecular dynamics

Release all position restraints, switch the barostat to the rigorous **Parrinello–Rahman** ensemble (`pcoupl = Parrinello-Rahman`, $\tau_p = 2.0\text{ ps}$), and launch production MD (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; source: Best Practices for Foundations in Molecular Simulations [Article.pdf):
```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
gmx mdrun -deffnm md_0_10
```
(source: gromacs_protein_ligand_complex_tutorial_lemkul.md)

## Step 7: Trajectory processing and quantitative analysis

### 1. PBC Re-wrapping and Structural Fitting (`trjconv`)
Raw MD trajectories display molecules broken across periodic boundaries. Correct this in two steps:
1. Re-wrap into compact molecular representations centered on the protein:
   ```bash
   gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
   ```
   Select `Protein` for centering and `System` for output (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).
2. Remove rigid-body rotation and translation by least-squares fitting onto the protein `Backbone`:
   ```bash
   gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o md_0_10_fit.xtc -fit rot+trans
   ```
   (source: gromacs_protein_ligand_complex_tutorial_lemkul.md)

### 2. Hydrogen-Bond Distance and Angle Geometry
Monitor specific non-covalent interactions (e.g., JZ4 phenolic hydroxyl `OAB`/`H12` to Gln102 sidechain `OE1`):
```bash
gmx distance -s md_0_10.tpr -f md_0_10_center.xtc -select 'resname "JZ4" and name OAB plus resid 102 and name OE1' -oall dist_hbond.xvg
```
Geometric criteria: donor–acceptor distance $\le 0.35\text{ nm}$ (observed: $0.31 \pm 0.05\text{ nm}$) and angle deviation $\le 30^\circ$ (observed: $23.3 \pm 9.4^\circ$), confirming a stable hydrogen bond (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

### 3. Ligand Heavy-Atom RMSD
Compute ligand RMSD relative to the energy-minimized starting structure:
```bash
gmx make_ndx -f em.gro -n index.ndx
# Select JZ4 heavy atoms -> name group "JZ4_Heavy"
gmx rms -s em.tpr -f md_0_10_center.xtc -tu ns -o rmsd_jz4.xvg
```
Select `Backbone` for least-squares fitting and `JZ4_Heavy` for RMSD calculation; a stable plateau at $\sim 0.1\text{ nm } (1.0\text{ \AA})$ indicates a stable binding pose and valid force-field parameterization (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

### 4. Nonbonded Interaction Energy Decomposition (`mdrun -rerun`)
Decompose short-range Coulombic (`Coul-SR`) and Lennard-Jones (`LJ-SR`) interaction energies between protein and ligand:
In `ie.mdp`, define:
```ini
energygrps = Protein JZ4
```
Re-evaluate energies across the trajectory using the CPU engine (`-nb cpu`, required by GROMACS when multiple `energygrps` are requested):
```bash
gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ie.tpr
gmx mdrun -deffnm ie -rerun md_0_10.xtc -nb cpu
gmx energy -f ie.edr -o interaction_energy.xvg
```
Observed CHARMM36 values: $\text{Coul-SR} = -20.5 \pm 7.4\text{ kJ/mol}$, $\text{LJ-SR} = -99.1 \pm 7.2\text{ kJ/mol}$, yielding a total gas-phase short-range interaction energy of $-119.6 \pm 10.3\text{ kJ/mol}$ (source: gromacs_protein_ligand_complex_tutorial_lemkul.md).

*(Thermodynamic Caveat: Gas-phase interaction energy is not a binding free energy $\Delta G_{\text{bind}}$; it omits ligand and receptor desolvation penalties as well as conformational and rotational/translational entropy losses. Quantitative affinity prediction requires [[mm-pbsa-and-mm-gbsa]], [[umbrella-sampling]], or [[alchemical-free-energy-methods]])* (source: gromacs_protein_ligand_complex_tutorial_lemkul.md; source: Best Practices for Foundations in Molecular Simulations [Article.pdf).

## Related pages
- [[lemkul-gromacs-protein-ligand-tutorial]]
- [[lemkul-2018-gromacs-tutorials]]
- [[molecular-dynamics-simulation]]
- [[force-fields]]
- [[open-force-field-and-ligand-parameterization]]
- [[simulation-and-benchmark-best-practices]]
- [[braun-2019-foundations-molecular-simulation]]
- [[ligand-preparation-and-protonation]]
- [[receptor-preparation]]
- [[binding-free-energy]]
- [[mm-pbsa-and-mm-gbsa]]

