"""Experiment planning calculations based on SAM001062 slide 25."""

from flex_crispr_probe_designer.constants import (
    CRISPR_REC_READS_PER_CELL,
    DEFAULT_CELLS_PER_GUIDE,
    EXPERIMENT_PLANNING,
    GEX_MIN_READS_PER_CELL,
)
from flex_crispr_probe_designer.models import ExperimentPlan


def plan_experiment(
    n_guides: int,
    cells_per_guide: int = DEFAULT_CELLS_PER_GUIDE,
) -> ExperimentPlan:
    """Calculate experiment requirements for a given guide library size.

    Uses the planning table from SAM001062 slide 25. For guide counts
    between table rows, rounds up to the next tier (conservative).
    For counts beyond the table, extrapolates from the largest tier.
    """
    min_cells = n_guides * cells_per_guide

    # Find the appropriate tier (round up)
    tier = EXPERIMENT_PLANNING[-1]  # Default to largest
    for row in EXPERIMENT_PLANNING:
        if n_guides <= row["guides"]:
            tier = row
            break

    # Scale values if beyond table range
    if n_guides > EXPERIMENT_PLANNING[-1]["guides"]:
        scale = n_guides / EXPERIMENT_PLANNING[-1]["guides"]
        wta_hybs = max(1, round(tier["wta_hybs"] * scale))
        gem_wells = max(1, round(tier["gem_wells"] * scale))
        cells_hybd = round(tier["cells_hybd_wta"] * scale)
        cells_fix = round(tier["cells_into_fixation"] * scale)
        cells_per_well = min_cells // gem_wells if gem_wells > 0 else min_cells
    else:
        wta_hybs = tier["wta_hybs"]
        gem_wells = tier["gem_wells"]
        cells_hybd = tier["cells_hybd_wta"]
        cells_fix = tier["cells_into_fixation"]
        cells_per_well = tier["cells_per_gem_well"]

    bc_oligo_hybs = min(tier["bc_oligo_hybs"], 48)

    gex_reads = min_cells * GEX_MIN_READS_PER_CELL
    crispr_reads = min_cells * CRISPR_REC_READS_PER_CELL

    return ExperimentPlan(
        n_guides=n_guides,
        cells_per_guide=cells_per_guide,
        min_cells=min_cells,
        wta_hybs=wta_hybs,
        bc_oligo_hybs=bc_oligo_hybs,
        gem_wells=gem_wells,
        cells_per_gem_well=cells_per_well,
        cells_hybd_wta=cells_hybd,
        cells_into_fixation=cells_fix,
        gex_reads_total=gex_reads,
        crispr_reads_total=crispr_reads,
        total_reads=gex_reads + crispr_reads,
    )
