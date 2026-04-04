"""Tests for Cell Ranger multi config generator."""

from pathlib import Path

import pytest

from flex_crispr_probe_designer.export.multiconfig import (
    build_multi_config,
    expand_well_range,
    generate_all_configs,
    parse_sample_sheet,
)
from flex_crispr_probe_designer.models import MultiConfigParams, SampleEntry


class TestExpandWellRange:
    def test_single_well(self) -> None:
        assert expand_well_range("A01") == ["A01"]

    def test_range_same_row(self) -> None:
        result = expand_well_range("A01:A04")
        assert result == ["A01", "A02", "A03", "A04"]

    def test_range_across_rows(self) -> None:
        result = expand_well_range("A11:B02")
        assert result == ["A11", "A12", "B01", "B02"]

    def test_comma_separated(self) -> None:
        result = expand_well_range("A01,A05,B03")
        assert result == ["A01", "A05", "B03"]

    def test_mixed(self) -> None:
        result = expand_well_range("A01:A03,B01")
        assert result == ["A01", "A02", "A03", "B01"]

    def test_full_row(self) -> None:
        result = expand_well_range("A01:A12")
        assert len(result) == 12
        assert result[0] == "A01"
        assert result[-1] == "A12"

    def test_invalid_well(self) -> None:
        with pytest.raises(ValueError, match="Invalid well"):
            expand_well_range("Z01")

    def test_invalid_column(self) -> None:
        with pytest.raises(ValueError, match="out of range"):
            expand_well_range("A13")

    def test_lowercase_normalized(self) -> None:
        assert expand_well_range("a01") == ["A01"]


class TestParseSampleSheet:
    def test_parse(self, tmp_path: Path) -> None:
        sheet = tmp_path / "samples.csv"
        sheet.write_text(
            "sample_id,plate,wells,gem_well,description\n"
            "WT,A,A01:A04,1,Wild type\n"
            "KO,A,A05:A08,1,Knockout\n"
        )
        entries = parse_sample_sheet(sheet)
        assert len(entries) == 2
        assert entries[0].sample_id == "WT"
        assert entries[0].wells == ["A01", "A02", "A03", "A04"]
        assert entries[0].gem_well == 1
        assert entries[1].sample_id == "KO"


class TestBuildMultiConfig:
    @pytest.fixture
    def params(self) -> MultiConfigParams:
        return MultiConfigParams(
            probe_set_path="/data/probe_set.csv",
            gex_fastqs_path="/data/fastqs/gex",
            crispr_fastqs_path="/data/fastqs/crispr",
            feature_ref_path="/data/feature_reference.csv",
        )

    @pytest.fixture
    def samples(self) -> list[SampleEntry]:
        return [
            SampleEntry(sample_id="WT", plate="A", wells=["A01", "A02"], gem_well=1),
            SampleEntry(sample_id="KO", plate="A", wells=["A03", "A04"], gem_well=1),
        ]

    def test_has_all_sections(self, samples: list[SampleEntry], params: MultiConfigParams) -> None:
        config = build_multi_config(samples, params)
        assert "[gene-expression]" in config
        assert "[libraries]" in config
        assert "[feature]" in config
        assert "[samples]" in config

    def test_probe_set_path(self, samples: list[SampleEntry], params: MultiConfigParams) -> None:
        config = build_multi_config(samples, params)
        assert "probe-set,/data/probe_set.csv" in config

    def test_create_bam_false(self, samples: list[SampleEntry], params: MultiConfigParams) -> None:
        config = build_multi_config(samples, params)
        assert "create-bam,false" in config

    def test_libraries(self, samples: list[SampleEntry], params: MultiConfigParams) -> None:
        config = build_multi_config(samples, params)
        assert "flex_gex,/data/fastqs/gex,Gene Expression" in config
        assert "flex_cr,/data/fastqs/crispr,CRISPR Guide Capture" in config

    def test_feature_ref(self, samples: list[SampleEntry], params: MultiConfigParams) -> None:
        config = build_multi_config(samples, params)
        assert "reference,/data/feature_reference.csv" in config

    def test_sample_barcodes(self, samples: list[SampleEntry], params: MultiConfigParams) -> None:
        config = build_multi_config(samples, params)
        assert "WT,A-A01|A-A02" in config
        assert "KO,A-A03|A-A04" in config


class TestGenerateAllConfigs:
    @pytest.fixture
    def params(self) -> MultiConfigParams:
        return MultiConfigParams(
            probe_set_path="/data/probe_set.csv",
            gex_fastqs_path="/data/fastqs/gex",
            crispr_fastqs_path="/data/fastqs/crispr",
            feature_ref_path="/data/feature_reference.csv",
        )

    def test_two_gem_wells(self, params: MultiConfigParams, tmp_path: Path) -> None:
        samples = [
            SampleEntry(sample_id="S1", plate="A", wells=["A01", "A02"], gem_well=1),
            SampleEntry(sample_id="S2", plate="A", wells=["A03", "A04"], gem_well=1),
            SampleEntry(sample_id="S3", plate="A", wells=["A05", "A06"], gem_well=2),
            SampleEntry(sample_id="S4", plate="A", wells=["A07", "A08"], gem_well=2),
        ]
        configs = generate_all_configs(samples, params, tmp_path / "out")
        assert len(configs) == 2
        assert (tmp_path / "out" / "multi_config_well1.csv").exists()
        assert (tmp_path / "out" / "multi_config_well2.csv").exists()

    def test_correct_sample_grouping(self, params: MultiConfigParams, tmp_path: Path) -> None:
        samples = [
            SampleEntry(sample_id="S1", plate="A", wells=["A01"], gem_well=1),
            SampleEntry(sample_id="S2", plate="A", wells=["A02"], gem_well=2),
        ]
        configs = generate_all_configs(samples, params, tmp_path / "out")
        # Well 1 should only contain S1
        assert "S1" in configs[0].config_text
        assert "S2" not in configs[0].config_text
        # Well 2 should only contain S2
        assert "S2" in configs[1].config_text
        assert "S1" not in configs[1].config_text


class TestPlateValidation:
    def test_valid_plates(self) -> None:
        for plate in ["A", "B", "C", "D"]:
            entry = SampleEntry(sample_id="x", plate=plate, wells=["A01"], gem_well=1)
            assert entry.plate == plate

    def test_invalid_plate(self) -> None:
        with pytest.raises(ValueError, match="Plate must be"):
            SampleEntry(sample_id="x", plate="E", wells=["A01"], gem_well=1)

    def test_lowercase_plate(self) -> None:
        entry = SampleEntry(sample_id="x", plate="a", wells=["A01"], gem_well=1)
        assert entry.plate == "A"
