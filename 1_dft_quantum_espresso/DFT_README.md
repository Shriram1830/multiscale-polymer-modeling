# DFT Modeling of the Crystalline PEK Unit Cell

This is the first stage of the DFT -> MD -> micromechanics pipeline. It relaxes the crystalline PEK unit cell and evaluates the resulting geometry, using the plane-wave pseudopotential package Quantum ESPRESSO. The output is used to evaluate the geometry adopted for the classical MD model in `2_md_lammps/md_crystalline/`, not to generate it directly - see Limitations below.

As with the rest of this repository, the files here are a vanilla template. Convergence ranges, cutoffs, and k-point grids reflect a reasonable testing procedure for a system of this size, but any specific numerical outcome (relaxed lattice parameters, convergence behavior, predicted density, etc.) belongs to a specific production run and is not reproduced here.

## Unit Cell Background

PEK adopts an orthorhombic unit cell in the space group Pbcn, with two monomers per cell (Z = 2), consistent with the broader PAEK family. The starting lattice parameters used in the template input files (`a = 7.62 Å`, `b = 5.87 Å`, `c = 10.03 Å`) are the experimentally reported values from Nishino, Tada, and Nakamae (*Polymer* 33, 1992, 736-743), used here as literature reference geometry, not as a result of this pipeline. The full unit cell contains 46 atoms (26 C, 16 H, 4 O).

Because no experimentally determined crystallographic information file (CIF) is available for PEK, `orthorhombic_cell_builder.py` generates approximate starting fractional coordinates programmatically: it takes a raw monomer, centers it, aligns the polymer backbone along the crystallographic c-axis, and applies the Pbcn symmetry operation to generate the second monomer, producing a CIF block suitable as a DFT starting structure.

## Workflow

1. **Build a starting structure** - `orthorhombic_cell_builder.py` produces placeholder fractional coordinates and a CIF block for an arbitrary monomer. Replace the dummy monomer coordinates with your own parameterized structure before use.

2. **SCF convergence testing** - `scf.in` is a template for testing plane-wave cutoff (`ecutwfc`) and k-point mesh convergence, with `_ECUTWFC_`, `_ECUTRHO_`, and `_KPOINTS_` placeholders meant to be swept by an external automation script. `degauss.in` is the equivalent template for testing Gaussian smearing width, with the `_DEGAUSS_` placeholder swept the same way.

3. **Variable-cell relaxation (`vc-relax`)** - once converged parameters are established, both the ionic positions and cell parameters are relaxed. This step is not included as a separate template file here since it reuses the same `&CONTROL`/`&SYSTEM` blocks as `scf.in` with `calculation = 'vc-relax'` and the relaxation settings described below.

## Convergence Testing Methodology

- **Plane-wave cutoff (`ecutwfc`)** - tested across a range of values; `scf.in` and `degauss.in` use 90 Ry as the illustrative converged value for this template.
- **Charge density cutoff (`ecutrho`)** - set to four times `ecutwfc`, as required for projector augmented-wave (PAW) pseudopotentials (360 Ry for the 90 Ry example above).
- **K-point mesh** - tested across several Monkhorst-Pack grids; a 6×6×3 mesh is used as the illustrative converged value.
- **Gaussian smearing width (`degauss`)** - tested across a range of values; 0.005 Ry is used as a conservatively small illustrative value in the absence of a clear convergence trend.

These are the values baked into the vanilla templates in this folder. For a real system, all four should be re-tested and the winning values substituted in place of the placeholders.

## Structural Relaxation Settings

The `vc-relax` step uses:

- **Exchange-correlation:** Perdew-Burke-Ernzerhof (PBE) generalized gradient approximation, with PAW pseudopotentials.
- **Dispersion correction:** Grimme's D3 semi-empirical dispersion correction (`vdw_corr = 'dft-d3'`) to capture long-range van der Waals interactions between chains.
- **Ionic and cell relaxation:** Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.
- **Cell degrees of freedom:** `cell_dofree = 'ibrav'` to relax only the a, b, and c lattice vector lengths while locking all cell angles at 90°, consistent with orthorhombic Pbcn symmetry.
- **Convergence thresholds:** thresholds on the SCF electronic loop, interatomic forces, total energy change between ionic steps, and cell pressure, each set to a value appropriate for a system of this size - see the `&CONTROL`/`&ELECTRONS` blocks in `scf.in` for the specific placeholder thresholds used in this template.

## Repository Contents

- `orthorhombic_cell_builder.py` - generates placeholder CIF for the crystalline unit cell
- `scf.in` - SCF template for cutoff/k-point convergence testing
- `degauss.in` - SCF template for smearing (degauss) convergence testing

## Tools Used

- **Quantum ESPRESSO** - plane-wave DFT relaxation and SCF convergence testing.
- **Python (NumPy)** - starting-structure generation in `orthorhombic_cell_builder.py`.

## Limitations and Assumptions

- No experimentally determined CIF is available for PEK, so the starting atomic coordinates are generated programmatically rather than taken from a crystallographic database. This is a known source of deviation from the true crystallographic energy minimum and should be treated as a starting approximation, not ground truth.
- The relaxed DFT geometry is used to evaluate the geometry adopted for the classical MD model in Stage 2, not to generate it directly. The MD unit cell in `2_md_lammps/md_crystalline/` is built independently from the literature lattice parameters and the programmatically generated coordinates.
- `ecutwfc`, `ecutrho`, k-point mesh, and `degauss` are all system-size- and pseudopotential-dependent and must be re-converged for any system other than the illustrative one used in this template.

## References

- Nishino, T., Tada, K. and Nakamae, K. *Polymer* 33 (1992): 736-743.
- Giannozzi, P. et al. *Journal of Physics: Condensed Matter* 21 (2009): 395502.
