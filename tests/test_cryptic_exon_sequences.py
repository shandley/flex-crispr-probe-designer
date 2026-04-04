"""Tests for pre-fetched cryptic exon probe sequences."""

from flex_crispr_probe_designer.design.cryptic_exon_sequences import (
    generate_all_cryptic_exon_probes,
    generate_stmn2_probes,
    generate_unc13a_probes,
)


class TestSTMN2Probes:
    def test_returns_three_probes(self) -> None:
        probes = generate_stmn2_probes()
        assert len(probes) == 3

    def test_probe_types(self) -> None:
        probes = generate_stmn2_probes()
        types = {p.probe_type for p in probes}
        assert types == {"inclusion_5prime", "inclusion_3prime", "skipping"}

    def test_junction_probe_length(self) -> None:
        probes = generate_stmn2_probes()
        for p in probes:
            assert len(p.junction_probe) == 50

    def test_gc_content_reasonable(self) -> None:
        probes = generate_stmn2_probes()
        for p in probes:
            assert 0.2 <= p.gc_content <= 0.8

    def test_inclusion_5prime_sequence(self) -> None:
        probes = generate_stmn2_probes()
        p5 = [p for p in probes if p.probe_type == "inclusion_5prime"][0]
        # Last 25 of upstream exon (CCTACAATGGCTAAAACAGCAATGG) + first 25 of CE
        assert "CAATGG" in p5.junction_probe  # Junction point from upstream exon
        assert "GACTCG" in p5.junction_probe  # Start of CE
        assert len(p5.junction_probe) == 50


class TestUNC13AProbes:
    def test_returns_three_probes(self) -> None:
        probes = generate_unc13a_probes()
        assert len(probes) == 3

    def test_all_50bp(self) -> None:
        probes = generate_unc13a_probes()
        for p in probes:
            assert len(p.junction_probe) == 50

    def test_gene_name(self) -> None:
        probes = generate_unc13a_probes()
        assert all(p.gene == "UNC13A" for p in probes)


class TestAllProbes:
    def test_generate_all(self) -> None:
        all_probes = generate_all_cryptic_exon_probes()
        assert "STMN2" in all_probes
        assert "UNC13A" in all_probes
        assert len(all_probes["STMN2"]) == 3
        assert len(all_probes["UNC13A"]) == 3
