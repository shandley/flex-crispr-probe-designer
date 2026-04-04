"""Tests for export modules."""

from pathlib import Path

import pandas as pd

from flex_crispr_probe_designer.design.probes import build_probe_set
from flex_crispr_probe_designer.export.cellranger import build_feature_reference, write_feature_reference
from flex_crispr_probe_designer.export.ordering import probes_to_opool_csv, probes_to_ordering_csv
from flex_crispr_probe_designer.export.spikein import calculate_spikein_pool, write_spikein_csv
from flex_crispr_probe_designer.models import GuideInput, OPoolScale, ScaffoldType


class TestOrderingCSV:
    def test_creates_file(self, sample_guides: list[GuideInput], tmp_path: Path) -> None:
        ps = build_probe_set(sample_guides, ScaffoldType.PREFERRED)
        out = probes_to_ordering_csv(ps, tmp_path / "ordering.csv")
        assert out.exists()
        df = pd.read_csv(out)
        # 1 LHS + 4 RHS = 5 rows
        assert len(df) == 5
        assert df.iloc[0]["Name"].startswith("LHS_")

    def test_rhs_has_phos(self, sample_guides: list[GuideInput], tmp_path: Path) -> None:
        ps = build_probe_set(sample_guides, ScaffoldType.PREFERRED)
        out = probes_to_ordering_csv(ps, tmp_path / "ordering.csv")
        df = pd.read_csv(out)
        for seq in df.iloc[1:]["Sequence"]:
            assert str(seq).startswith("/5Phos/")


class TestOPoolCSV:
    def test_creates_file(self, sample_guides: list[GuideInput], tmp_path: Path) -> None:
        ps = build_probe_set(sample_guides, ScaffoldType.PREFERRED)
        out = probes_to_opool_csv(ps, tmp_path / "opool.csv")
        assert out.exists()
        df = pd.read_csv(out)
        assert len(df) == 4  # Only RHS probes
        assert "Pool Name" in df.columns


class TestFeatureReference:
    def test_build_rows(self, sample_guides: list[GuideInput]) -> None:
        rows = build_feature_reference(sample_guides)
        assert len(rows) == 4
        assert all(r.feature_type == "CRISPR Guide Capture" for r in rows)
        assert all(r.read == "R2" for r in rows)

    def test_pattern_preferred(self, sample_guides: list[GuideInput]) -> None:
        rows = build_feature_reference(sample_guides, ScaffoldType.PREFERRED)
        # Verified against 10X's own 1.2M K562 Flex CRISPR feature reference
        assert rows[0].pattern == "GCTATGCTGTTTCCAGCTTAGCTCTTAAAC(BC)"

    def test_pattern_at_scaffold(self, sample_guides: list[GuideInput]) -> None:
        rows = build_feature_reference(sample_guides, ScaffoldType.AT_SUBSTITUTION)
        assert rows[0].pattern == "CTTGCTATGCTGTTTCCAGCATAGCTCTTAAAC(BC)"

    def test_non_targeting(self, sample_guides: list[GuideInput]) -> None:
        rows = build_feature_reference(sample_guides)
        nt = [r for r in rows if r.target_gene_id == "Non-Targeting"]
        assert len(nt) == 1
        assert nt[0].target_gene_name == "Non-Targeting"

    def test_writes_csv(self, sample_guides: list[GuideInput], tmp_path: Path) -> None:
        out = write_feature_reference(sample_guides, tmp_path / "feature_ref.csv")
        assert out.exists()
        df = pd.read_csv(out)
        assert list(df.columns) == [
            "id", "name", "read", "pattern", "sequence",
            "feature_type", "target_gene_id", "target_gene_name",
        ]


class TestSpikeIn:
    def test_preferred_scaffold(self) -> None:
        pool = calculate_spikein_pool(ScaffoldType.PREFERRED)
        assert pool.lhs_final_conc_nm == 400
        assert pool.lhs_vol_ul == 2.2
        assert pool.total_vol_ul == 44.0

    def test_at_scaffold_increases_lhs(self) -> None:
        pool = calculate_spikein_pool(ScaffoldType.AT_SUBSTITUTION)
        assert pool.lhs_final_conc_nm == 1200
        assert pool.lhs_vol_ul == 6.6
        # IDTE reduced by 4.4 µL
        assert pool.idte_vol_ul < 39.6
        assert pool.total_vol_ul == 44.0

    def test_10pmol_scale(self) -> None:
        pool = calculate_spikein_pool(ScaffoldType.PREFERRED, OPoolScale.PMOL_10)
        assert pool.rhs_vol_ul == 11.0

    def test_writes_csv(self, tmp_path: Path) -> None:
        pool = calculate_spikein_pool(ScaffoldType.PREFERRED)
        out = write_spikein_csv(pool, tmp_path / "spikein.csv")
        assert out.exists()
