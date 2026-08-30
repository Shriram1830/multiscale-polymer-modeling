"""
crystalline_tensile_properties.py

Parses LAMMPS log files from uniaxial tensile deformation simulations of 
crystalline structures. Extracts axial stress/strain and lateral strains to 
simultaneously calculate Young's Modulus and Poisson's Ratio. 
Generates combined subplots and a final directional (orthotropic) property table.
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
STRAIN_MIN = 0.01
STRAIN_MAX = 0.03

# Column indices in your LAMMPS log (0-indexed).
STRAIN_COLS = {"X": 13, "Y": 14, "Z": 15}
STRESS_COLS  = {"X": 19, "Y": 20, "Z": 21}

TRANSVERSE = {
    "X": ("Y", "Z"),
    "Y": ("X", "Z"),
    "Z": ("X", "Y"),
}

DIRECTIONS = ["X", "Y", "Z"]

def read_log(filepath, axial_dir):
    t1_dir, t2_dir = TRANSVERSE[axial_dir]
    axial_strain, axial_stress, t1_strain, t2_strain = [], [], [], []
    reading = False

    # Encoding fix included to prevent UnicodeDecodeError crashes
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("Loop time"):
                reading = False
            if line.startswith("Per MPI rank"):
                reading = True
                continue
            if reading:
                cols = line.split()
                if len(cols) >= 22:
                    try:
                        axial_strain.append(float(cols[STRAIN_COLS[axial_dir]]))
                        axial_stress.append(float(cols[STRESS_COLS[axial_dir]]))
                        t1_strain.append(float(cols[STRAIN_COLS[t1_dir]]))
                        t2_strain.append(float(cols[STRAIN_COLS[t2_dir]]))
                    except ValueError:
                        pass

    return (np.array(axial_strain), np.array(axial_stress),
            np.array(t1_strain),    np.array(t2_strain))

print("=" * 70)
print("PROCESSING CRYSTALLINE TENSION LOGS (MODULUS & POISSON'S RATIO)")
print("=" * 70 + "\n")

summary_data = []

# Dictionary to store E and nu separated by axis for anisotropic averaging
axis_data = {
    "X": {"E": [], "nu": []},
    "Y": {"E": [], "nu": []},
    "Z": {"E": [], "nu": []}
}

for seed in range(1, N_SEEDS + 1):
    for direction in DIRECTIONS:

        # Standardized vanilla filename
        log_file = f"crystalline_seed{seed}_tension_{direction}.log"

        if not os.path.exists(log_file):
            print(f"[SKIP] {log_file} not found")
            continue

        eps, sig, t1, t2 = read_log(log_file, direction)

        if len(eps) == 0:
            print(f"[WARN] {log_file} — no data rows parsed")
            continue

        mask = (eps >= STRAIN_MIN) & (eps <= STRAIN_MAX)
        if mask.sum() < 3:
            print(f"[WARN] {log_file} — fewer than 3 points in elastic window")
            continue

        # -------------------------------------------------------------
        # 1. Young's Modulus Calculation (Linear Fit)
        # -------------------------------------------------------------
        slope_E, intercept_E = np.polyfit(eps[mask], sig[mask], 1)
        E_gpa = slope_E / 1000.0

        # -------------------------------------------------------------
        # 2. Poisson's Ratio Calculation (Linear Fit)
        # -------------------------------------------------------------
        slope_nu1, intercept_nu1 = np.polyfit(eps[mask], t1[mask], 1)
        slope_nu2, intercept_nu2 = np.polyfit(eps[mask], t2[mask], 1)

        nu1 = -slope_nu1
        nu2 = -slope_nu2
        nu_avg = (nu1 + nu2) / 2.0

        t1_name, t2_name = TRANSVERSE[direction]

        # Store properties for the combined table
        summary_data.append({
            "Seed": seed,
            "Axis": direction,
            "E (GPa)": round(E_gpa, 3),
            f"ν ({direction}{t1_name})": round(nu1, 4),
            f"ν ({direction}{t2_name})": round(nu2, 4),
            "Average ν": round(nu_avg, 4)
        })

        # Store properties per axis for final directional averaging
        axis_data[direction]["E"].append(E_gpa)
        axis_data[direction]["nu"].append(nu_avg)

        # =========================================================
        # COMBINED PLOTTING (1x2 Subplots)
        # =========================================================
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(f"Crystalline Polymer Tensile Test - Seed {seed} {direction}", fontsize=14, fontweight='bold')

        # --- Subplot 1: Stress-Strain (Young's Modulus) ---
        ax1.plot(eps, sig, label=f"Stress-Strain", color='black')
        fit_eps = np.array([STRAIN_MIN, STRAIN_MAX])
        fit_sig_E = slope_E * fit_eps + intercept_E
        ax1.plot(fit_eps, fit_sig_E, 'r--', linewidth=2, label=f"Fit (E = {E_gpa:.2f} GPa)")
        
        ax1.set_xlabel("Engineering Strain")
        ax1.set_ylabel("Stress (MPa)")
        ax1.set_title("Young's Modulus")
        ax1.legend()
        ax1.grid(alpha=0.3)

        # --- Subplot 2: Poisson's Ratio (Lateral vs Axial Strain) ---
        ax2.scatter(eps[mask], t1[mask], alpha=0.5, s=25, label=f"{direction} vs {t1_name}", color='steelblue', marker='o')
        ax2.scatter(eps[mask], t2[mask], alpha=0.5, s=25, label=f"{direction} vs {t2_name}", color='darkseagreen', marker='s')

        ax2.plot(fit_eps, slope_nu1 * fit_eps + intercept_nu1, color='blue', linewidth=2, label=f"Fit (ν = {nu1:.4f})")
        ax2.plot(fit_eps, slope_nu2 * fit_eps + intercept_nu2, color='green', linewidth=2, label=f"Fit (ν = {nu2:.4f})")

        ax2.set_xlabel(f"Axial Strain (v_etrue{direction.lower()})")
        ax2.set_ylabel("Lateral Strain")
        ax2.set_title(f"Poisson's Ratio (Avg ν = {nu_avg:.4f})")
        ax2.legend()
        ax2.grid(alpha=0.3)

        # Save and close to prevent memory warnings
        plt.tight_layout()
        save_name = f"plot_crystal_seed{seed}_mechanics_{direction}.png"
        plt.savefig(save_name, dpi=150)
        plt.close()

# =========================================================
# FINAL OUTPUT TABLE AND DIRECTIONAL SUMMARY
# =========================================================
print("\n" + "=" * 70)
print("    CRYSTALLINE POLYMER — COMBINED RESULTS TABLE")
print("=" * 70)

df_summary = pd.DataFrame(summary_data)
display(df_summary)

print("\n" + "=" * 70)
print("    DIRECTIONAL AVERAGES (PER AXIS)")
print("=" * 70)

# Calculate and print stats per axis (Critical for Anisotropic Materials)
for direction in DIRECTIONS:
    E_list = axis_data[direction]["E"]
    nu_list = axis_data[direction]["nu"]

    if len(E_list) > 0:
        E_arr = np.array(E_list)
        nu_arr = np.array(nu_list)

        E_mean = E_arr.mean()
        E_std = E_arr.std()
        nu_mean = nu_arr.mean()
        nu_std = nu_arr.std()

        print(f"--- AXIS {direction} ---")
        print(f"  E (Young's modulus) = {E_mean:.3f} ± {E_std:.3f} GPa")
        print(f"  ν (Poisson's ratio) = {nu_mean:.3f} ± {nu_std:.3f}\n")
    else:
        print(f"--- AXIS {direction} ---")
        print("  No data available.\n")
