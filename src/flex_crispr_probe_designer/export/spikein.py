"""Spike-in pool preparation calculations."""

from pathlib import Path

import pandas as pd

from flex_crispr_probe_designer.constants import (
    LHS_CONC_AT_SUBST_NM,
    LHS_CONC_PREFERRED_NM,
    LHS_VOL_AT_SCAFFOLD_UL,
    RHS_CONC_PER_PROBE_NM,
    SPIKEIN_VOLUMES,
)
from flex_crispr_probe_designer.models import OPoolScale, ScaffoldType, SpikeInPool


def calculate_spikein_pool(
    scaffold: ScaffoldType,
    opool_scale: OPoolScale = OPoolScale.PMOL_50,
    n_hybridizations: int = 16,
) -> SpikeInPool:
    """Calculate spike-in pool preparation volumes.

    For A>T scaffold, LHS volume is increased to 6.6 µL (for 1200 nM)
    and IDTE is reduced proportionally.
    """
    base_volumes = SPIKEIN_VOLUMES[opool_scale.value]

    lhs_vol = base_volumes["lhs_8um"]
    idte_vol = base_volumes["idte_ph8"]
    rhs_vol = base_volumes["rhs_stock"]
    total = base_volumes["total"]

    if scaffold == ScaffoldType.AT_SUBSTITUTION:
        # Increase LHS from 2.2 to 6.6 µL, reduce IDTE by the difference
        delta = LHS_VOL_AT_SCAFFOLD_UL - lhs_vol
        lhs_vol = LHS_VOL_AT_SCAFFOLD_UL
        idte_vol = idte_vol - delta
        lhs_conc = LHS_CONC_AT_SUBST_NM
    else:
        lhs_conc = LHS_CONC_PREFERRED_NM

    return SpikeInPool(
        opool_scale=opool_scale,
        scaffold_type=scaffold,
        n_hybridizations=n_hybridizations,
        idte_vol_ul=round(idte_vol, 1),
        rhs_vol_ul=round(rhs_vol, 1),
        lhs_vol_ul=round(lhs_vol, 1),
        total_vol_ul=round(total, 1),
        lhs_final_conc_nm=lhs_conc,
        rhs_final_conc_nm=RHS_CONC_PER_PROBE_NM,
    )


def write_spikein_csv(pool: SpikeInPool, output: Path) -> Path:
    """Write spike-in pool preparation table as CSV."""
    rows = [
        {"Component": "IDTE, pH 8.0", "Volume (µL)": pool.idte_vol_ul},
        {"Component": "RHS CRISPR, resuspended stock", "Volume (µL)": pool.rhs_vol_ul},
        {"Component": f"LHS CRISPR, 8 µM ({pool.lhs_final_conc_nm} nM final)", "Volume (µL)": pool.lhs_vol_ul},
        {"Component": "Total", "Volume (µL)": pool.total_vol_ul},
    ]
    df = pd.DataFrame(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    return output
