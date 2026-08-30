"""
crystalline_shear_properties.py

Parses LAMMPS log files from shear deformation simulations.
Extracts shear stress and strain, calculates the Shear Modulus (G) 
via linear regression in the elastic regime, and generates plots.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from IPython.display import display

# =====================================================================
# CONFIGURATION
# =====================================================================
N_SEEDS = 5
STRAIN_MIN = 0.04
STRAIN_MAX = 0.06

# Column indices for SHEAR in your LAMMPS log (0-indexed).
STRAIN_COLS = {"XY": 11, "XZ": 12, "YZ": 13} # v_eengxy, v_eengxz, v_eengyz
STRESS_COLS = {"XY": 27, "XZ": 28, "YZ": 29} # f_sxy_ave, f_sxz_ave, f_syz_ave

DIRECTIONS = ["XY", "YZ", "XZ"]

def read_shear_log(filepath, shear_dir):
    shear_strain, shear_stress = [], []
    reading = False

    # Using errors="ignore" to prevent the UnicodeDecodeError!
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("Loop time"):
                reading = False
            if line.startswith("Per MPI rank"):
                reading = True
                continue
            if reading:
                cols = line.split()
                if len(cols) >= 30:
                    try:
                        # Extract engineering shear strain
                        shear_strain.append(float(cols[STRAIN_COLS[shear_dir]]))

                        # Extract stress and convert from atm to MPa (1 atm = 0.101325 MPa)
                        stress_atm = float(cols[STRESS_COLS[shear_dir]])
                        shear_stress.append(stress_atm * 0.101325)
                    except ValueError:
                        pass

    return np.array(shear_strain), np.array(shear_stress)


print("=" * 60)
print("PROCESSING INDIVIDUAL SHEAR FILES")
print("=" * 60 + "\n")

summary_data = []
all_G = []

for seed in range(1, N_SEEDS + 1):
    for direction in DIRECTIONS:

        # Standardized file naming convention
        log_file = f"crystalline_seed{seed}_shear_{direction}.log"

        if not os.path.exists(log_file):
            print(f"[SKIP] {log_file} not found")
            continue

        eps, sig = read_shear_log(log_file, direction)

        if len(eps) == 0:
            print(f"[WARN] {log_file} — no data rows parsed")
            continue

        mask = (eps >= STRAIN_MIN) & (eps <= STRAIN_MAX)
        if mask.sum() < 3:
            print(f"[WARN] {log_file} — fewer than 3 points in elastic window")
            continue

        # Calculate Shear Modulus (G) from the linear fit of Shear Stress vs Strain
        slope, intercept = np.polyfit(eps[mask], sig[mask], 1)
        G_gpa = slope / 1000.0  # MPa to GPa

        # Store properties for the summary table
        summary_data.append({
            "Seed": seed,
            "Shear Plane": direction,
            "Shear Modulus G (GPa)": round(G_gpa, 3)
        })
        all_G.append(G_gpa)

        # =========================================================
        # GENERATE PLOT
        # =========================================================
        plt.figure(figsize=(8, 6))

        # Plot Raw Data
        plt.plot(eps, sig, color='blue', label=f"Shear Stress-Strain ({direction})")

        # Plot Fit Line
        fit_eps = np.array([STRAIN_MIN, STRAIN_MAX])
        fit_sig = slope * fit_eps + intercept
        plt.plot(fit_eps, fit_sig, 'r--', linewidth=2, label=f"Elastic Fit (G = {G_gpa:.2f} GPa)")

        plt.xlabel("Engineering Shear Strain")
        plt.ylabel("Shear Stress (MPa)")
        plt.title(f"Polymer Shear Test | Seed {seed} | Plane {direction}")
        plt.legend()
        plt.grid(alpha=0.3)

        # Save image automatically
        save_name = f"plot_seed{seed}_shear_{direction}.png"
        plt.savefig(save_name, dpi=150, bbox_inches="tight")
        plt.close()

# =========================================================
# FINAL OUTPUT TABLE
# =========================================================
print("\n" + "=" * 60)
print("    POLYMER — SHEAR MODULUS TABLE")
print("=" * 60)

if summary_data:
    df_summary = pd.DataFrame(summary_data)
    display(df_summary)
    
    # Calculate and print overall statistics
    G_arr = np.array(all_G)
    print("\n" + "-" * 60)
    print("  OVERALL ISOTROPIC SHEAR MODULUS (All Planes & Seeds)")
    print("-" * 60)
    print(f"  G (Shear modulus) = {G_arr.mean():.3f} ± {G_arr.std():.3f} GPa")
    print("=" * 60)
else:
    print("No valid data processed. Please ensure your files are uploaded correctly.")
