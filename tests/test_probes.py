"""Tests for probe generation."""

import pytest

from flex_crispr_probe_designer.constants import (
    AT_SCAFFOLD_RC_DNA,
    AT_SCAFFOLD_RC_IDT,
    PREFERRED_SCAFFOLD_RC,
    RHS_LINKER,
    RHS_PCONSTANT_SEQ,
    TRUSEQ_R2_HANDLE,
)
from flex_crispr_probe_designer.design.dna import reverse_complement
from flex_crispr_probe_designer.design.probes import build_lhs_probe, build_probe_set, build_rhs_probe
from flex_crispr_probe_designer.models import CustomScaffold, GuideInput, ScaffoldType


class TestBuildLHSProbe:
    def test_preferred_scaffold(self) -> None:
        lhs = build_lhs_probe(ScaffoldType.PREFERRED)
        assert lhs.sequence == TRUSEQ_R2_HANDLE + PREFERRED_SCAFFOLD_RC
        assert lhs.length == len(TRUSEQ_R2_HANDLE) + len(PREFERRED_SCAFFOLD_RC)
        assert lhs.has_lna is False
        assert lhs.idt_notation == lhs.sequence  # No mods for preferred

    def test_at_scaffold(self) -> None:
        lhs = build_lhs_probe(ScaffoldType.AT_SUBSTITUTION)
        assert lhs.sequence == TRUSEQ_R2_HANDLE + AT_SCAFFOLD_RC_DNA
        assert lhs.has_lna is True
        assert AT_SCAFFOLD_RC_IDT in lhs.idt_notation
        assert "+" in lhs.idt_notation

    def test_preferred_length(self) -> None:
        lhs = build_lhs_probe(ScaffoldType.PREFERRED)
        assert lhs.length == 34 + 25  # TruSeq R2 (34) + scaffold RC (25)

    def test_at_length(self) -> None:
        lhs = build_lhs_probe(ScaffoldType.AT_SUBSTITUTION)
        assert lhs.length == 34 + 28  # TruSeq R2 (34) + AT scaffold RC (28)


class TestBuildRHSProbe:
    def test_structure(self, single_guide: GuideInput) -> None:
        rhs = build_rhs_probe(single_guide)
        # Structure: TAAAC + spacer_RC + CCCATATAAGAAA
        assert rhs.sequence.startswith(RHS_LINKER)
        assert rhs.sequence.endswith(RHS_PCONSTANT_SEQ)
        assert rhs.length == len(RHS_LINKER) + 20 + len(RHS_PCONSTANT_SEQ)

    def test_idt_notation_has_phos(self, single_guide: GuideInput) -> None:
        rhs = build_rhs_probe(single_guide)
        assert rhs.idt_notation.startswith("/5Phos/")

    def test_spacer_rc(self, single_guide: GuideInput) -> None:
        rhs = build_rhs_probe(single_guide)
        # GAATGGGAACGTACGGGAAA -> RC = TTTCCCGTACGTTCCCATTC
        assert rhs.spacer_rc == "TTTCCCGTACGTTCCCATTC"
        assert rhs.spacer_rc in rhs.sequence

    def test_rhs_length(self, single_guide: GuideInput) -> None:
        rhs = build_rhs_probe(single_guide)
        assert rhs.length == 38  # 5 + 20 + 13


class TestBuildProbeSet:
    def test_probe_set(self, sample_guides: list[GuideInput]) -> None:
        ps = build_probe_set(sample_guides, ScaffoldType.PREFERRED)
        assert ps.n_guides == 4
        assert len(ps.rhs_probes) == 4
        assert ps.lhs.scaffold_type == ScaffoldType.PREFERRED

    def test_all_rhs_unique(self, sample_guides: list[GuideInput]) -> None:
        ps = build_probe_set(sample_guides, ScaffoldType.PREFERRED)
        sequences = [r.sequence for r in ps.rhs_probes]
        assert len(sequences) == len(set(sequences))


class TestCustomScaffold:
    def test_custom_scaffold_basic(self) -> None:
        # Use a 25 bp scaffold (same length as preferred)
        scaffold_seq = "AGAGCTAAGCTGGAAACAGCATAGC"  # The preferred scaffold itself
        custom = CustomScaffold(scaffold_sequence=scaffold_seq)
        lhs = build_lhs_probe(ScaffoldType.CUSTOM, custom)

        expected_rc = reverse_complement(scaffold_seq)
        assert lhs.scaffold_rc == expected_rc
        assert lhs.sequence.endswith(expected_rc)
        assert lhs.scaffold_type == ScaffoldType.CUSTOM
        assert lhs.has_lna is False
        assert "custom" in lhs.scaffold_label

    def test_custom_scaffold_with_idt_notation(self) -> None:
        scaffold_seq = "AGAGCTAAGCTGGAAACAGCATAGC"
        idt = "+GCTATGCTGTTTC+CAGCTTAGCTCT"
        custom = CustomScaffold(scaffold_sequence=scaffold_seq, idt_notation=idt)
        lhs = build_lhs_probe(ScaffoldType.CUSTOM, custom)

        assert lhs.has_lna is True
        assert idt in lhs.idt_notation

    def test_custom_scaffold_different_length(self) -> None:
        # A shorter scaffold (e.g., 20 bp)
        scaffold_seq = "AGAGCTAAGCTGGAAACAGC"
        custom = CustomScaffold(scaffold_sequence=scaffold_seq)
        lhs = build_lhs_probe(ScaffoldType.CUSTOM, custom)

        assert lhs.length == 34 + 20  # TruSeq R2 + 20 bp scaffold RC
        assert "20 bp" in lhs.scaffold_label

    def test_custom_scaffold_in_probe_set(self, sample_guides: list[GuideInput]) -> None:
        scaffold_seq = "AGAGCTAAGCTGGAAACAGCATAGC"
        custom = CustomScaffold(scaffold_sequence=scaffold_seq)
        ps = build_probe_set(sample_guides, ScaffoldType.CUSTOM, custom)

        assert ps.scaffold_type == ScaffoldType.CUSTOM
        assert ps.lhs.scaffold_type == ScaffoldType.CUSTOM
        assert ps.n_guides == 4

    def test_custom_requires_scaffold_object(self) -> None:
        with pytest.raises(ValueError, match="custom scaffold must be provided"):
            build_lhs_probe(ScaffoldType.CUSTOM)

    def test_custom_scaffold_validates_dna(self) -> None:
        with pytest.raises(ValueError, match="non-ACGT"):
            CustomScaffold(scaffold_sequence="ACGTNNN")
