# Multiscale Modeling of Poly Ether Ketone (PEK)

This repository holds the modeling pipeline for predicting the bulk mechanical properties of semicrystalline Polyetherketone (PEK) across three length scales: quantum-mechanical, atomistic, and microstructural. It accompanies the research paper *"Multiscale Modeling of Poly Ether Ketone (PEK) using DFT, MD and Micromechanics."*

All code in this repository is a **vanilla template**: input files use placeholder geometry, lattice parameters, and force-field constants, and analysis scripts use illustrative default values. Nothing here reproduces the specific numerical results reported in the paper. The intent is to document a reusable, general-purpose workflow, not to publish the underlying dataset.

## Motivation

Polyetherketone is a semicrystalline thermoplastic from the polyaryletherketone (PAEK) family, used in place of metals in high-temperature aerospace structural components. Like other PAEK resins, its crystalline and amorphous regions coexist across a hierarchy of length scales - from the unit cell up through granular crystal blocks, lamella stacks, and spherulites - and the degree of crystallinity is a process-dependent variable rather than a fixed material constant.

No single simulation method spans this entire range: molecular dynamics (MD) can predict single-phase elastic constants at the nanometer scale, but the spherulitic microstructure that governs bulk behavior exists at the micron scale. A multiscale framework is needed to connect the two.

This repository implements that framework: density functional theory (DFT) relaxes and evaluates the crystalline unit cell, MD predicts the elastic constants of the crystalline and amorphous phases independently, and a hierarchical micromechanics homogenization assembles those phase-level constants into bulk semicrystalline properties.

## Modeling Framework

The pipeline has three modeling stages, run in order of increasing length scale, followed by a shared data-analysis layer:

1. **DFT (Quantum ESPRESSO)** - relaxes the crystalline PEK unit cell starting from literature lattice parameters and programmatically generated fractional coordinates, and evaluates the resulting geometry against the experimentally reported crystal structure.

2. **Molecular Dynamics (LAMMPS)** - predicts elastic constants for the crystalline and amorphous phases independently, each modeled as its own MD system with its own equilibration, tensile, and (for the crystalline phase) shear deformation protocol.

3. **Micromechanics (ANSYS Material Designer)** - homogenizes the two phase-level stiffness tensors up through three explicit length scales (Granular Crystal Block -> Lamella -> Spherulite) to predict bulk semicrystalline properties as a function of crystallinity.

4. **Data Analysis (Python)** - parses LAMMPS log output from stage 2 into density, modulus, and Poisson's ratio values and generates the plots and summary tables consumed by stage 3.

## Repository Structure

    .
    ├── 1_dft_quantum_espresso/
    │   ├── orthorhombic_cell_builder.py  # generates placeholder CIF for the crystalline unit cell
    │   ├── scf.in                        # SCF template for cutoff/k-point convergence testing
    │   └── degauss.in                    # SCF template for smearing (degauss) convergence testing
    │
    ├── 2_md_lammps/
    │   ├── md_crystalline/
    │   │   ├── equilibration.in           # energy minimization + NVT/NPT equilibration
    │   │   ├── tensile.in                 # uniaxial tensile deformation, x/y/z
    │   │   ├── triclinic.in               # orthogonal -> triclinic box conversion
    │   │   └── shear.in                   # shear deformation, xy/xz/yz
    │   │
    │   └── md_amorphous/
    │       └── equilibration.in           # energy minimization + NVT/NPT equilibration
    │
    ├── 3_micromechanics_fea_ansys/
    │   ├── MICROMECHANICS_README.md       # stage-specific documentation
    │   ├── lamella_length_calculator.py   # crystallinity -> amorphous plank length
    │   └── spherulite_builder.py          # NxNxN spherulite RUC generator
    │
    └── 4_data_analysis/
        ├── density_calculation.py             # equilibrium density from LAMMPS logs
        ├── crystalline_tensile_properties.py  # E and nu from crystalline tensile logs
        ├── crystalline_shear_properties.py   # G from crystalline shear logs
        └── amorphous_properties.py            # isotropic E and nu from amorphous tensile logs

## Stage 1 - DFT (Quantum ESPRESSO)

Structural relaxation and convergence testing on the crystalline PEK unit cell. The plane-wave cutoff, charge-density cutoff, k-point mesh, and Gaussian smearing width are all convergence-tested before the production variable-cell relaxation (`vc-relax`) is run.

The relaxed geometry is evaluated against the literature crystal structure, not carried forward automatically - the classical MD model in Stage 2 is built independently from the literature lattice parameters and programmatically generated coordinates, using the tools listed below.

See `1_dft_quantum_espresso/DFT_README.md` for the full methodology.

## Stage 2 - Molecular Dynamics (LAMMPS)

The crystalline and amorphous phases are modeled as two independent systems:

- **Crystalline:** a single orthorhombic unit cell (not replicated into a supercell, to avoid breaking bonding connectivity across periodic images) is built externally in Avogadro, VESTA, and OVITO, then parameterized with the PCFF-IFF force field via the LUNAR toolkit before being handed to LAMMPS for equilibration, uniaxial tensile testing along all three axes, and shear testing in all three principal planes (via an orthogonal-to-triclinic box conversion).

- **Amorphous:** built via in-situ reactive polymerization (LAMMPS `fix bond/react`) from randomly packed monomer templates, followed by densification, high-temperature annealing, and equilibration, then characterized by uniaxial tensile testing along all three axes and averaged into a single isotropic response.

Multiple independent replicate models (different initial velocity seeds) are used for both phases to support statistical averaging.

See `2_md_lammps/MD_README.md` for the full methodology.

## Stage 3 - Micromechanics (ANSYS Material Designer)

Homogenizes the two MD-derived phase stiffness tensors up through the GCB -> Lamella -> Spherulite hierarchy.

See `3_micromechanics_fea_ansys/MICROMECHANICS_README.md` for the full methodology, the two supporting Python scripts, and their illustrative default values.

## Stage 4 - Data Analysis

Generic LAMMPS log parsers that extract equilibrium density (`density_calculation.py`), tensile modulus and Poisson's ratio (`crystalline_tensile_properties.py`, `amorphous_properties.py`), and shear modulus (`crystalline_shear_properties.py`) via linear regression over a specified elastic-strain window, and generate the corresponding stress-strain plots.

## Naming Convention

Files are named `{phase}_seed{n}_{test}_{axis}.{ext}`, where `phase` is `crystalline` or `amorphous`.

Each MD input script's output filename is what the corresponding data-analysis script expects as input:

- `crystalline_seed1_equilibrated.data`
- `crystalline_seed1_triclinic.data`
- `crystalline_seed1_tension_X.log` -> `crystalline_tensile_properties.py`
- `crystalline_seed1_shear_XY.log` -> `crystalline_shear_properties.py`
- `amorphous_seed1_equilibrated.data`
- `amorphous_seed1_tension_X.log` -> `amorphous_properties.py`

## Tools and Dependencies

- **Quantum ESPRESSO** - DFT relaxation and convergence testing
- **Avogadro, VESTA, OVITO** - structure building and inspection
- **LUNAR toolkit** - force field atom typing and LAMMPS data file generation
- **LAMMPS** - all molecular dynamics simulation
- **ANSYS Material Designer** - all micromechanics homogenization
- **Python** (NumPy, pandas, Matplotlib) - geometry calculations and log-file analysis

## Status and Limitations

- All input files in this repository use placeholder values (lattice parameters, force-field constants, target crystallinities, etc.) and are meant to be adapted to a specific system, not run as-is.
- The DFT and MD stages do not yet have their own per-folder READMEs; this document and the inline comments in each script are the current documentation for those stages.
- See `3_micromechanics_fea_ansys/MICROMECHANICS_README.md` for stage-specific limitations and assumptions in the micromechanics homogenization.

## References

1. Mallick, P.K. *Fiber-Reinforced Composites: Materials, Manufacturing, and Design.* CRC Press, 2007.
2. Wang, Y. et al. *RSC Advances* 6 (2016): 3198-3209.
3. Talbott, M.F., Springer, G.S., Berglund, L.A. *Journal of Composite Materials* 21 (1987): 1056-1081.
4. Bandyopadhyay, A. et al. *Polymer* 52 (2011): 2445-2452.
5. Pisani, W. A. et al. *Polymer* 163 (2019): 96-105.
6. Kashmari, K. et al. *ACS Applied Engineering Materials* 1 (2023): 3167-3177.
7. Kemppainen, J. et al. *Journal of Chemical Information and Modeling* 64 (2024): 5108-5126.
8. Blundell, D.J. and Osborn, B.N. *Polymer* 24 (1983): 953-958.
9. Nishino, T., Tada, K. and Nakamae, K. *Polymer* 33 (1992): 736-743.
10. Heinz, H. et al. *Langmuir* 29 (2013): 1754-1765.
11. Thompson, A.P. et al. *Computer Physics Communications* 271 (2022): 108171.
12. Martyna, G.J., Tobias, D.J. and Klein, M.L. *Journal of Chemical Physics* 101 (1994): 4177-4189.
13. Giannozzi, P. et al. *Journal of Physics: Condensed Matter* 21 (2009): 395502.
14. ANSYS, Inc. *ANSYS Material Designer.* ANSYS, Inc., Canonsburg, PA, 2023.
