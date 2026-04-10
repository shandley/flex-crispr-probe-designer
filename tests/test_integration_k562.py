"""Integration test: validate probe designer output against real 10X K562 Flex CRISPR data."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from flex_crispr_probe_designer.export.cellranger import build_feature_reference
from flex_crispr_probe_designer.models import GuideInput, ScaffoldType

K562_FEATURE_REF = Path(
    "/Users/scotthandley/Code/tools/10X_flex_v2/flex_v2_analysis"
    "/data/10x_flex_crispr_k562/count_feature_reference.csv"
)

pytestmark = pytest.mark.skipif(
    not K562_FEATURE_REF.exists(),
    reason=f"K562 feature reference not found at {K562_FEATURE_REF}",
)

DNA_20BP = re.compile(r"^[ACGT]{20}$")


@pytest.fixture(scope="module")
def k562_crispr_rows() -> pd.DataFrame:
    """Load only CRISPR Guide Capture rows from the K562 feature reference."""
    df = pd.read_csv(K562_FEATURE_REF)
    crispr = df[df["feature_type"] == "CRISPR Guide Capture"].copy()
    assert len(crispr) > 0, "No CRISPR Guide Capture rows found in K562 feature reference"
    return crispr


class TestK562PatternMatch:
    """Verify the probe designer produces the same pattern as the real K562 dataset."""

    def test_pattern_matches_k562(self, k562_crispr_rows: pd.DataFrame) -> None:
        """The PREFERRED scaffold pattern must match all K562 CRISPR Guide Capture rows."""
        k562_patterns = k562_crispr_rows["pattern"].unique()
        assert len(k562_patterns) == 1, (
            f"Expected a single pattern across all K562 CRISPR rows, got {len(k562_patterns)}"
        )
        k562_pattern = k562_patterns[0]

        # Build a feature reference with a dummy guide using PREFERRED scaffold
        dummy = [GuideInput(guide_id="test", spacer_sequence="AAAACCCCGGGGTTTTAAAA", target_gene="TEST")]
        rows = build_feature_reference(dummy, ScaffoldType.PREFERRED)
        designer_pattern = rows[0].pattern

        assert designer_pattern == k562_pattern, (
            f"Pattern mismatch!\n"
            f"  Designer:  {designer_pattern}\n"
            f"  K562 data: {k562_pattern}"
        )

    def test_k562_feature_type(self, k562_crispr_rows: pd.DataFrame) -> None:
        """All CRISPR rows must have feature_type 'CRISPR Guide Capture'."""
        assert (k562_crispr_rows["feature_type"] == "CRISPR Guide Capture").all()

    def test_k562_sequences_are_valid_20bp_dna(self, k562_crispr_rows: pd.DataFrame) -> None:
        """Every CRISPR Guide Capture sequence must be a valid 20-bp DNA string."""
        for idx, seq in k562_crispr_rows["sequence"].items():
            assert DNA_20BP.match(seq), (
                f"Row {idx}: sequence '{seq}' is not a valid 20-bp DNA string"
            )

    def test_k562_read_field(self, k562_crispr_rows: pd.DataFrame) -> None:
        """All CRISPR Guide Capture rows should specify R2."""
        assert (k562_crispr_rows["read"] == "R2").all()
