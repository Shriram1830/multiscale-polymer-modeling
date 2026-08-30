"""
spherulite_builder.py

Generates the subcell-level crystallinity assignment for a 3D cubic 
Repeating Unit Cell (RUC) representing a semicrystalline polymer spherulite.
Outputs a CSV for import into finite element homogenization software 
(e.g., ANSYS Material Designer).

Algorithm summary:
1. Divides the RUC grid into concentric shells based on squared radial distance.
2. Optimizes the assignment of discrete crystallinity levels to these shells to 
   match a target overall crystallinity.
3. Enforces a strict monotonically decreasing crystallinity gradient from core to edge.
4. Calculates the radial unit vector (chain-axis orientation) for each subcell.
"""

from dataclasses import dataclass
from itertools import product
import csv

# Discrete lamella crystallinities available for assignment
CRYSTALLINITY_LEVELS_PCT: list[float] = [75, 60, 45, 30, 15, 0]

GRID_SIZE: int = 6                  # 6x6x6 = 216 subcells
TARGET_OVERALL_CRYSTALLINITY_PCT: float = 30.0


@dataclass
class Subcell:
    """One subcell of the spherulite RUC."""
    x: int                  # grid coordinate, 1-indexed (1..GRID_SIZE)
    y: int
    z: int
    r_squared: float        # squared distance from cube center
    crystallinity_pct: float
    radial_dir: tuple[float, float, float]  # unit vector, chain-axis orientation


def _cube_center(grid_size: int) -> float:
    """Center coordinate for a 1-indexed grid of the given size."""
    return (grid_size + 1) / 2.0


def _build_r2_shells(grid_size: int) -> tuple[list[float], dict[float, int]]:
    """
    Groups subcells into shells based on squared radial distance to preserve symmetry.
    """
    center = _cube_center(grid_size)
    r2_counts: dict[float, int] = {}

    for x, y, z in product(range(1, grid_size + 1), repeat=3):
        dx, dy, dz = x - center, y - center, z - center
        r2 = dx**2 + dy**2 + dz**2
        r2_counts[r2] = r2_counts.get(r2, 0) + 1

    r2_sorted = sorted(r2_counts.keys())
    return r2_sorted, r2_counts


def _search_closest_symmetric_assignment(
    r2_sorted: list[float],
    r2_counts: dict[float, int],
    levels_pct: list[float],
    target_pct: float,
    total_cells: int,
    require_all_levels: bool = True,
) -> list[int]:
    """
    Exhaustively searches monotonically non-increasing assignments of crystallinity
    levels to radial shells to find the volume-weighted average closest to the target.
    """
    n_shells = len(r2_sorted)
    n_levels = len(levels_pct)
    counts = [r2_counts[r2] for r2 in r2_sorted]
    target_sum = target_pct * total_cells

    best_seq: list[int] | None = None
    best_diff = float("inf")

    def recurse(shell_idx: int, min_level_idx: int, seq: list[int]) -> None:
        nonlocal best_seq, best_diff

        if shell_idx == n_shells:
            if require_all_levels and len(set(seq)) < n_levels:
                return
            total = sum(counts[s] * levels_pct[seq[s]] for s in range(n_shells))
            diff = abs(total - target_sum)
            if diff < best_diff:
                best_diff = diff
                best_seq = list(seq)
            return

        if require_all_levels:
            remaining_shells = n_shells - shell_idx
            remaining_levels_needed = n_levels - 1 - min_level_idx
            if remaining_shells < remaining_levels_needed:
                return

        for level_idx in range(min_level_idx, n_levels):
            seq.append(level_idx)
            recurse(shell_idx + 1, level_idx, seq)
            seq.pop()

    recurse(0, 0, [])

    if best_seq is None:
        raise RuntimeError(
            "No valid shell assignment found. Try require_all_levels=False."
        )
    return best_seq


def build_spherulite(
    grid_size: int = GRID_SIZE,
    levels_pct: list[float] = None,
    target_pct: float = TARGET_OVERALL_CRYSTALLINITY_PCT,
    require_all_levels: bool = True,
) -> list[Subcell]:
    """
    Builds the full list of subcells with assigned crystallinity and radial orientation.
    """
    if levels_pct is None:
        levels_pct = CRYSTALLINITY_LEVELS_PCT

    center = _cube_center(grid_size)
    total_cells = grid_size ** 3

    r2_sorted, r2_counts = _build_r2_shells(grid_size)
    shell_assignment = _search_closest_symmetric_assignment(
        r2_sorted, r2_counts, levels_pct, target_pct, total_cells,
        require_all_levels=require_all_levels,
    )
    r2_to_level: dict[float, float] = {
        r2: levels_pct[shell_assignment[i]] for i, r2 in enumerate(r2_sorted)
    }

    subcells: list[Subcell] = []
    for x, y, z in product(range(1, grid_size + 1), repeat=3):
        dx, dy, dz = x - center, y - center, z - center
        r2 = dx**2 + dy**2 + dz**2

        norm = (dx**2 + dy**2 + dz**2) ** 0.5
        radial_dir = (dx / norm, dy / norm, dz / norm) if norm > 1e-9 else (1.0, 0.0, 0.0)

        subcells.append(
            Subcell(
                x=x, y=y, z=z,
                r_squared=r2,
                crystallinity_pct=r2_to_level[r2],
                radial_dir=radial_dir,
            )
        )

    return subcells


def overall_crystallinity_pct(subcells: list[Subcell]) -> float:
    """Volume-weighted average crystallinity across all subcells."""
    return sum(c.crystallinity_pct for c in subcells) / len(subcells)


def write_csv(subcells: list[Subcell], filepath: str) -> None:
    """Exports the subcell table to CSV format."""
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Subcube_X", "Subcube_Y", "Subcube_Z", "r_squared",
            "Crystallinity_pct", "Radial_dir_x", "Radial_dir_y", "Radial_dir_z",
        ])
        for c in subcells:
            writer.writerow([
                c.x, c.y, c.z, c.r_squared, c.crystallinity_pct,
                round(c.radial_dir[0], 4),
                round(c.radial_dir[1], 4),
                round(c.radial_dir[2], 4),
            ])


if __name__ == "__main__":
    subcells = build_spherulite()

    avg = overall_crystallinity_pct(subcells)
    print(f"Total subcells: {len(subcells)}")
    print(f"Overall average crystallinity: {avg:.4f}%  "
          f"(target: {TARGET_OVERALL_CRYSTALLINITY_PCT}%)")

    print("\nCell count per crystallinity level:")
    counts_by_level: dict[float, int] = {}
    for c in subcells:
        counts_by_level[c.crystallinity_pct] = counts_by_level.get(c.crystallinity_pct, 0) + 1
    for level in sorted(counts_by_level, reverse=True):
        print(f"  {level:5.1f}%  ->  {counts_by_level[level]:3d} subcells")

    output_file = "Polymer_Spherulite_6x6x6.csv"
    write_csv(subcells, output_file)
    print(f"\nWrote {output_file}")