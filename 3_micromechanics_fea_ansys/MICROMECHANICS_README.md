# Micromechanics Modeling of Semi-Crystalline PEK

This part of the project predicts the bulk mechanical properties of
semi-crystalline Polyetherketone (PEK) by homogenizing molecular-scale
elastic properties up through three hierarchical length scales —
Granular Crystal Block (GCB) → Lamella → Spherulite — using ANSYS
Material Designer. It is the final stage of a larger DFT → MD →
micromechanics pipeline; this document covers only the micromechanics
stage.

## Why Micromechanics Is Needed

Molecular dynamics gives you the elastic properties of two idealized,
single-phase materials: 100% crystalline PEK and 100% amorphous PEK.
Neither one is what you'd actually measure in a lab. Real PEK is a
semi-crystalline solid where ordered crystalline lamellae and
disordered amorphous regions coexist at the microstructural level, and
the ratio between them depends on processing history. Getting from
"two idealized phases" to "one real bulk material" requires a
homogenization step that accounts for how those two phases are
geometrically arranged and how load is shared between them. That's
what this stage does.

We followed the same overall framework as Pisani et al. (2019),
*"Multiscale modeling of PEEK using reactive molecular dynamics
modeling and micromechanics,"* Polymer 163, 96–105, adapting it from
PEEK to PEK and from their MAC/GMC-based approach to ANSYS Material
Designer.

## Hierarchical Structure

Semi-crystalline PAEK-family polymers (PEK included) organize into
four nested structural levels, from largest to smallest: root network,
spherulite, lamella stack, and granular crystal block. Following
Pisani et al., we did not model the root network explicitly — a single
spherulite surrounded by amorphous material is treated as
representative of the bulk semi-crystalline microstructure.

### Level 1 — Granular Crystal Block (GCB)

The GCB is the smallest structural unit: primary crystalline blocks
bound together by a less-ordered "secondary structure," which is
modeled as amorphous PEK. Following the GCB dimensions reported by
Wang et al. (*RSC Advances* 6, 2016, 3198–3209), the GCB was built as
a 120 nm cube with an internal crystalline volume fraction of 86.2%.

This was modeled in ANSYS Material Designer as a two-phase RVE
combining the MD-derived crystalline PEK stiffness tensor
(orthotropic) with the MD-derived amorphous PEK stiffness (isotropic),
using a predefined short-fiber-composite geometry template at the
correct volume fraction. The output is a single effective, still
mildly anisotropic, GCB stiffness tensor — this becomes the input
material for Level 2.

### Level 2 — Lamella

A lamella is modeled as a plank: one GCB block placed end-to-end with
a variable-length amorphous PEK block along the chain axis. Because
the GCB length is fixed (120 nm), the *overall* crystallinity of a
given plank is entirely controlled by how long the attached amorphous
block is. Longer amorphous section → lower plank crystallinity.

The relationship comes from two nested volume fractions. The
crystallinity of a plank isn't the same as the volume fraction of GCB
material in that plank, because the GCB itself is only 86.2%
crystalline internally:

```
V_crys  = V_GCB × 0.862
V_GCB   = L_GCB / (L_GCB + L_amorphous)
```

Solving for the amorphous length given a target crystallinity:

```
L_amorphous = L_GCB × (1 / V_GCB − 1),   where V_GCB = V_crys / 0.862
```

We built five lamella planks at 15%, 30%, 45%, 60%, and 75%
crystallinity, each requiring a specific amorphous block length:

| Target crystallinity | GCB volume fraction | Amorphous length |
|---|---|---|
| 15% | 0.174 | 569.6 nm |
| 30% | 0.348 | 224.8 nm |
| 45% | 0.522 | 109.9 nm |
| 60% | 0.696 | 52.4 nm |
| 75% | 0.870 | 17.9 nm |

Each plank was built and solved as its own two-body RVE in ANSYS
Material Designer with periodic boundary conditions, GCB material from
Level 1 assigned to one block and amorphous PEK to the other, giving
five separate effective lamella stiffness tensors. Note that ANSYS's
embedded geometry kernel enforces a minimum modeling scale in the
micrometer range, so in practice all geometry was built at a
proportionally scaled-up size (nanometers → millimeters, a 10⁶×
scale factor) rather than at the true physical scale. Linear elastic
homogenization is scale-invariant — the effective stiffness tensor
depends only on volume fraction and shape ratio, not absolute
dimensions — so this has no effect on the resulting properties. This
is a common and defensible workaround for CAD-kernel precision limits
in RVE-based micromechanics, not a physical approximation.

`lamella_length_calculator.py` in this repository performs the
crystallinity-to-length calculation above for an arbitrary target list.

### Level 3 — Spherulite

The spherulite is where the model becomes genuinely three-dimensional.
Physically, lamellae in a spherulite radiate outward from a central
nucleus, and crystallinity is highest at the core (where lamellae pack
densely) and decreases toward the outer edge (where they've spread
apart and amorphous material dominates).

We captured this with a 6×6×6 cubic RUC (216 subcells). Each subcell
was assigned:

1. **A crystallinity level**, drawn from the five lamella materials
   built in Level 2 (plus pure amorphous PEK for the outermost
   subcells), based on the subcell's squared radial distance from the
   cube center. Crystallinity is strictly non-increasing as radial
   distance increases, and every subcell at the same radial distance
   receives the identical crystallinity level, preserving full
   spherical symmetry.

2. **A radial orientation vector**, pointing from the cube center to
   that subcell's center. This defines the local chain-axis direction
   for the anisotropic lamella material in that subcell, following
   the same logic used by Pisani et al. (their Fig. 10): each subcell's
   lamella "points" outward from the spherulite nucleus, so the
   assembled material has no single preferred direction at the bulk
   scale — a defining feature of a real spherulite.

Because ANSYS Material Designer works with a small, discrete set of
pre-solved lamella materials rather than a continuously graded
property field, matching a specific target crystallinity exactly isn't
always possible — the achievable overall averages form a fixed,
irregular set of numbers depending on how many subcells fall into each
radial shell. We treated this as a constrained search problem: among
every radially-symmetric, monotonically-decreasing assignment of the
five crystallinity levels (plus 0%) to the grid's nine distinct radial
shells, find the one whose volume-weighted average crystallinity comes
closest to the typical experimentally observed bulk value of 30%
(Hudson et al., *Macromolecules* 25, 1992, 1759–1765), while requiring
that all five lamella materials actually appear at least once.

The resulting spherulite averages **29.44% crystallinity**, distributed
as follows:

| Crystallinity | Subcell count |
|---|---|
| 75% | 8 |
| 60% | 24 |
| 45% | 56 |
| 30% | 48 |
| 15% | 24 |
| 0% (amorphous) | 56 |
| **Total** | **216** |

`spherulite_builder.py` generates this assignment from scratch — it
enumerates the radial shell structure of a 6×6×6 grid, searches for the
best crystallinity-to-shell mapping under the symmetry and monotonicity
constraints described above, and exports a CSV with each subcell's
coordinates, crystallinity, and radial orientation vector, ready to
drive the corresponding ANSYS setup.

## Tools Used

- **ANSYS Material Designer** — all RVE homogenization at every level
  (GCB, lamella, spherulite), using periodic boundary conditions
  throughout.
- **Python** — geometry/crystallinity calculations, the discrete shell
  optimization search, and CSV export for the spherulite subcell table.

## Repository Contents

```
lamella_length_calculator.py   # crystallinity -> amorphous plank length
spherulite_builder.py          # 6x6x6 spherulite RUC generator
PEK_Spherulite_6x6x6.csv       # generated 216-subcell output
```

## Limitations and Assumptions

- The GCB internal crystalline volume fraction (86.2%) and root-network
  simplification (single spherulite representing the bulk) are carried
  over from the PEEK literature (Wang et al., Pisani et al.) rather
  than measured independently for PEK, since PEK and PEEK share the
  same PAEK-family GCB morphology.
- The spherulite is modeled as a cube rather than a sphere, and the
  lamella as a rectangular plank rather than the rounded shapes seen
  in SEM imaging — both are standard simplifications in RUC-based
  micromechanics that preserve the correct volume fractions and
  overall topology while remaining computationally tractable.
- Interfacial effects between crystalline and amorphous regions are
  not explicitly modeled at the molecular level, consistent with
  Pisani et al.'s treatment of PEEK.
- Geometry was built at a proportionally scaled-up size to work around
  ANSYS's CAD kernel minimum feature size; this does not affect the
  resulting homogenized properties (see Level 2 above).

## References

- Pisani, W.A. et al. "Multiscale modeling of PEEK using reactive
  molecular dynamics modeling and micromechanics." *Polymer* 163
  (2019): 96–105.
- Wang, Y. et al. "Unusual crystalline morphology of poly aryl ether
  ketones (PAEKs)." *RSC Advances* 6 (2016): 3198–3209.
- Hudson, S.D. et al. "Semicrystalline morphology of poly(aryl ether
  ether ketone)/poly(ether imide) blends." *Macromolecules* 25 (1992):
  1759–1765.
