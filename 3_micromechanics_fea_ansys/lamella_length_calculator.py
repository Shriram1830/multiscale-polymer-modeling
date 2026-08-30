"""
lamella_length_calculator.py
Computes the amorphous plank length required to achieve a target overall 
crystallinity for a two-block lamella RUC in ANSYS Material Designer.
"""

from dataclasses import dataclass

# Material parameters (Placeholders for vanilla template)
GCB_INTERNAL_CRYSTAL_FRACTION: float = 0.850
DEFAULT_GCB_LENGTH_NM: float = 100.0


@dataclass
class LamellaPlank:
    target_crystallinity_pct: float
    gcb_length_nm: float
    gcb_volume_fraction: float
    amorphous_length_nm: float


def compute_amorphous_length(
    target_crystallinity_pct: float,
    gcb_length_nm: float = DEFAULT_GCB_LENGTH_NM,
    gcb_internal_fraction: float = GCB_INTERNAL_CRYSTAL_FRACTION,
) -> LamellaPlank:
    """
    Calculates the amorphous block length needed to hit a target crystallinity.
    """
    v_crys = target_crystallinity_pct / 100.0
    v_gcb = v_crys / gcb_internal_fraction

    if v_gcb > 1.0:
        max_pct = gcb_internal_fraction * 100.0
        raise ValueError(
            f"Target {target_crystallinity_pct}% is impossible. "
            f"Max for this GCB is {max_pct:.2f}%."
        )

    amorphous_length_nm = gcb_length_nm * (1.0 / v_gcb - 1.0)

    return LamellaPlank(
        target_crystallinity_pct=target_crystallinity_pct,
        gcb_length_nm=gcb_length_nm,
        gcb_volume_fraction=v_gcb,
        amorphous_length_nm=amorphous_length_nm,
    )


def build_lamella_table(target_crystallinities_pct: list[float]) -> list[LamellaPlank]:
    """Generates plank geometries for a list of target crystallinities."""
    return [compute_amorphous_length(pct) for pct in target_crystallinities_pct]


if __name__ == "__main__":
    targets_pct = [20, 30, 40, 50, 60]
    planks = build_lamella_table(targets_pct)

    print(f"{'Target %':<12}{'GCB Vf':<12}{'L_amorphous (nm)':<20}")
    print("-" * 44)
    for p in planks:
        print(
            f"{p.target_crystallinity_pct:<12.1f}"
            f"{p.gcb_volume_fraction:<12.5f}"
            f"{p.amorphous_length_nm:<20.2f}"
        )