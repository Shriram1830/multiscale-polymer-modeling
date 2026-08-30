"""
density_calculation.py

Calculates the equilibrium density from LAMMPS log files by averaging the 
last N steps (default: 1000) of the thermodynamic output.
Usage: python density_calculation *.log
"""

import sys, math, glob, os

LAST_N = 1000
DENSITY_COLS = ["density"]

def mean(vals): 
    return sum(vals) / len(vals)

def stdev(vals):
    if len(vals) < 2: return 0.0
    m = mean(vals)
    return math.sqrt(sum((x - m)**2 for x in vals) / (len(vals) - 1))

def parse_log(filepath):
    density_vals = []
    col_idx = None
    reading = False

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.split()
            if not parts: continue
            
            lower_parts = [p.lower() for p in parts]
            if any(name in lower_parts for name in DENSITY_COLS):
                col_idx = next(i for i, c in enumerate(lower_parts) if c in DENSITY_COLS)
                reading = True
                continue

            if reading and col_idx is not None:
                try:
                    int(parts[0]) # Validates step number
                    density_vals.append(float(parts[col_idx]))
                except (ValueError, IndexError):
                    reading = False

    return density_vals

def main():
    args = sys.argv[1:]
    if not args:
        print("Usage: python lammps_density_calc.py *.log")
        sys.exit(1)

    filepaths = sorted([f for arg in args for f in glob.glob(arg)])
    if not filepaths:
        print("No files found.")
        sys.exit(1)

    print(f"\nDensity Calc (Last {LAST_N} steps)")
    print("-" * 65)
    print(f"{'File':<30} {'N':>6} {'Mean (g/cc)':>12} {'Std':>12}")
    print("-" * 65)

    seed_means = []
    for fpath in filepaths:
        if not os.path.isfile(fpath): continue
        
        vals = parse_log(fpath)[-LAST_N:]
        if not vals:
            print(f"{os.path.basename(fpath):<30} No data")
            continue

        m, s = mean(vals), stdev(vals)
        seed_means.append(m)
        fname = os.path.basename(fpath)
        print(f"{fname[:28]:<30} {len(vals):>6} {m:>12.5f} {s:>12.5f}")

    if seed_means:
        print("-" * 65)
        print(f"Overall: {mean(seed_means):.4f} ± {stdev(seed_means):.4f} g/cc ({len(seed_means)} files)\n")

if __name__ == "__main__":
    main()
