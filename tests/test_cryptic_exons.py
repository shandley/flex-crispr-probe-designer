"""Tests for cryptic exon probe design module."""

from flex_crispr_probe_designer.design.cryptic_exons import (
    CRYPTIC_EXON_CATALOG,
    design_cryptic_exon_probes_from_sequences,
    design_junction_probe,
    get_cryptic_exon,
    list_cryptic_exons,
)


class TestCatalog:
    def test_catalog_not_empty(self) -> None:
        assert len(CRYPTIC_EXON_CATALOG) >= 5

    def test_list_returns_copy(self) -> None:
        a = list_cryptic_exons()
        b = list_cryptic_exons()
        assert a is not b

    def test_get_stmn2(self) -> None:
        ce = get_cryptic_exon("STMN2")
        assert ce is not None
        assert ce.gene == "STMN2"
        assert ce.chrom == "chr8"
        assert ce.confidence == "HIGH"
        assert ce.ce_length == 227

    def test_get_unc13a(self) -> None:
        ce = get_cryptic_exon("UNC13A")
        assert ce is not None
        assert ce.strand == "-"
        assert ce.ce_length == 128

    def test_get_case_insensitive(self) -> None:
        assert get_cryptic_exon("stmn2") is not None
        assert get_cryptic_exon("Stmn2") is not None

    def test_get_nonexistent(self) -> None:
        assert get_cryptic_exon("FAKE_GENE") is None

    def test_all_have_required_fields(self) -> None:
        for ce in CRYPTIC_EXON_CATALOG:
            assert ce.gene
            assert ce.chrom.startswith("chr")
            assert ce.ce_start < ce.ce_end
            assert ce.strand in ("+", "-")
            assert ce.confidence in ("HIGH", "MEDIUM", "LOW")


class TestDesignJunctionProbe:
    def test_basic_probe(self) -> None:
        probe = design_junction_probe(
            ce_sequence_at_junction="A" * 25,
            flanking_sequence_at_junction="C" * 25,
            gene="TEST",
            probe_type="inclusion",
        )
        assert len(probe.junction_probe) == 50
        assert probe.junction_probe == "C" * 25 + "A" * 25
        assert probe.gene == "TEST"
        assert probe.probe_type == "inclusion"

    def test_gc_content(self) -> None:
        probe = design_junction_probe(
            ce_sequence_at_junction="ACGTACGTACGTACGTACGTACGTG",
            flanking_sequence_at_junction="ACGTACGTACGTACGTACGTACGTG",
            gene="TEST",
        )
        assert 0.4 <= probe.gc_content <= 0.6

    def test_rc_is_correct(self) -> None:
        probe = design_junction_probe(
            ce_sequence_at_junction="A" * 25,
            flanking_sequence_at_junction="C" * 25,
            gene="TEST",
        )
        assert probe.junction_probe_rc == "T" * 25 + "G" * 25


class TestDesignFromSequences:
    def test_returns_three_probes(self) -> None:
        probes = design_cryptic_exon_probes_from_sequences(
            gene="TEST",
            upstream_exon_3prime="A" * 30,
            ce_5prime="C" * 30,
            ce_3prime="G" * 30,
            downstream_exon_5prime="T" * 30,
        )
        assert len(probes) == 3
        types = {p.probe_type for p in probes}
        assert types == {"inclusion_5prime", "inclusion_3prime", "skipping"}

    def test_inclusion_5prime(self) -> None:
        probes = design_cryptic_exon_probes_from_sequences(
            gene="TEST",
            upstream_exon_3prime="A" * 30,
            ce_5prime="C" * 30,
            ce_3prime="G" * 30,
            downstream_exon_5prime="T" * 30,
        )
        p5 = [p for p in probes if p.probe_type == "inclusion_5prime"][0]
        # Should be last 25 of upstream + first 25 of CE
        assert p5.junction_probe == "A" * 25 + "C" * 25

    def test_skipping_probe(self) -> None:
        probes = design_cryptic_exon_probes_from_sequences(
            gene="TEST",
            upstream_exon_3prime="A" * 30,
            ce_5prime="C" * 30,
            ce_3prime="G" * 30,
            downstream_exon_5prime="T" * 30,
        )
        skip = [p for p in probes if p.probe_type == "skipping"][0]
        # Should be last 25 of upstream + first 25 of downstream (no CE)
        assert skip.junction_probe == "A" * 25 + "T" * 25
