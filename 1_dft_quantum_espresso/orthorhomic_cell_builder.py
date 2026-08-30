import numpy as np

# =============================================================================
# Orthorhombic Polymer Crystal Coordinate Generator
# Space group: Pbcn (No. 60) | Z = 2
#
# Description: 
#   Takes a raw Cartesian monomer, centers it, aligns the polymer backbone 
#   along the crystallographic c-axis, and applies the Pbcn symmetry operation 
#   (2_1 screw axis) to generate the second monomer in the unit cell.
#   Outputs a formatted CIF file block.
# =============================================================================

# ── 1. Lattice Parameters (Placeholder values) ────────────────────────────────
a = 10.00   # Å  (orthorhombic a-axis)
b = 8.00    # Å  (orthorhombic b-axis)
c = 12.00   # Å  (chain axis — each monomer spans approximately c/2 = 6.00 Å)

# ── 2. Raw Monomer Coordinates ────────────────────────────────────────────────
# INPUT: Replace with your parameterized monomer coordinates.
# Define the indices that represent the repeat vector of the polymer chain.
TAIL_INDEX = 0  # Starting atom of the repeating backbone segment
HEAD_INDEX = 4  # Ending atom of the repeating backbone segment

# Dummy 5-atom polymer chain for demonstration
monomer_raw = [
    ("C", 1.0, 1.0, 0.0),  # idx 0 (Tail)
    ("C", 1.5, 1.2, 1.5),  # idx 1
    ("O", 1.0, 1.0, 3.0),  # idx 2
    ("C", 0.5, 0.8, 4.5),  # idx 3
    ("C", 1.0, 1.0, 6.0),  # idx 4 (Head)
]

elements = [row[0] for row in monomer_raw]
coords   = np.array([row[1:] for row in monomer_raw])

# ── 3. Center at origin ───────────────────────────────────────────────────────
coords_centered = coords - np.mean(coords, axis=0)

# ── 4. Align the chain axis to crystallographic c (z) ────────────────────────
# The bridge vector from TAIL -> HEAD is the polymer repeat direction. 
# Rotating this onto z ensures the chain runs along c.
bridge_vec = coords_centered[HEAD_INDEX] - coords_centered[TAIL_INDEX]
bridge_len = np.linalg.norm(bridge_vec)
bridge_vec /= bridge_len

print(f"Bridge vector length (Tail → Head): {bridge_len:.4f} Å")
print(f"Half unit cell c/2:                 {c/2:.4f} Å")
print(f"Difference (chain strain):          {bridge_len - c/2:.4f} Å\n")

z_axis = np.array([0.0, 0.0, 1.0])
v      = np.cross(bridge_vec, z_axis)
c_dot  = np.dot(bridge_vec, z_axis)
s_val  = np.linalg.norm(v)

# Rotation matrix using Rodrigues' rotation formula
if s_val > 1e-6:
    vx = np.array([[ 0,     -v[2],  v[1]],
                   [ v[2],   0,    -v[0]],
                   [-v[1],   v[0],  0   ]])
    R = np.eye(3) + vx + vx @ vx * ((1 - c_dot) / s_val**2)
else:
    R = np.eye(3) if c_dot > 0 else -np.eye(3)

aligned = coords_centered @ R.T   # (N, 3) — chain now along z

# ── 5. Translate monomer 1 to (a/4, b/4, c/4) in Cartesian space ─────────────
aligned[:, 0] += a * 0.25
aligned[:, 1] += b * 0.25
aligned[:, 2] += c * 0.25

# ── 6. Fractional coordinates of monomer 1 ───────────────────────────────────
fx1 = (aligned[:, 0] / a) % 1.0
fy1 = (aligned[:, 1] / b) % 1.0
fz1 = (aligned[:, 2] / c) % 1.0

# ── 7. Generate monomer 2 via Pbcn symmetry operation ─────────────────────────
# Pbcn general position #2:  (x, y, z) → (−x, −y, z + ½)
# This represents a 180-degree rotation about the c-axis and a c/2 translation.
fx2 = (-fx1) % 1.0         # −x mod 1
fy2 = (-fy1) % 1.0         # −y mod 1
fz2 = (fz1 + 0.5) % 1.0    #  z + ½

# ── 8. Sanity check: Min distance between monomer 1 and monomer 2 ─────────────
print("Sanity check — Minimum intermolecular distance:")
min_dist = np.inf
for i in range(len(fx1)):
    for j in range(len(fx2)):
        # Apply periodic boundary conditions
        dx = min(abs(fx2[j] - fx1[i]), 1 - abs(fx2[j] - fx1[i])) * a
        dy = min(abs(fy2[j] - fy1[i]), 1 - abs(fy2[j] - fy1[i])) * b
        dz = min(abs(fz2[j] - fz1[i]), 1 - abs(fz2[j] - fz1[i])) * c
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        if dist < min_dist:
            min_dist = dist
print(f"  Closest atomic distance: {min_dist:.3f} Å\n")

# ── 9. Print the CIF atom_site loop & full file ───────────────────────────────
print("=" * 60)
print("COMPLETE CIF FILE — copy everything below into polymer_crystal.cif")
print("=" * 60)
print(f"""
data_polymer_crystal

_cell_length_a    {a:.4f}
_cell_length_b    {b:.4f}
_cell_length_c    {c:.4f}
_cell_angle_alpha  90.000
_cell_angle_beta   90.000
_cell_angle_gamma  90.000

_symmetry_space_group_name_H-M  'P b c n'
_symmetry_Int_Tables_number      60

loop_
_symmetry_equiv_pos_as_xyz
'x, y, z'
'-x, -y, z+1/2'
'-x, y+1/2, -z+1/2'
'x, -y+1/2, -z'
'-x, -y, -z'
'x, y, -z+1/2'
'x, -y+1/2, z+1/2'
'-x, y+1/2, z'

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z""")

counter = 1
for el, x, y, z in zip(elements, fx1, fy1, fz1):
    print(f"{el}{counter}  {x:.6f}  {y:.6f}  {z:.6f}")
    counter += 1
for el, x, y, z in zip(elements, fx2, fy2, fz2):
    print(f"{el}{counter}  {x:.6f}  {y:.6f}  {z:.6f}")
    counter += 1