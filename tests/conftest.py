"""Shared test fixtures."""

from __future__ import annotations

import pytest

from flex_crispr_probe_designer.models import GuideInput


@pytest.fixture
def sample_guides() -> list[GuideInput]:
    """A small guide library for testing."""
    return [
        GuideInput(guide_id="sg_BRCA1_1", spacer_sequence="GAATGGGAACGTACGGGAAA", target_gene="BRCA1"),
        GuideInput(guide_id="sg_BRCA1_2", spacer_sequence="TCCTCAAGGAGCTTTATTGA", target_gene="BRCA1"),
        GuideInput(guide_id="sg_TP53_1", spacer_sequence="CCATTGTTCAATATCGTCCG", target_gene="TP53"),
        GuideInput(guide_id="sg_NT_1", spacer_sequence="ACGGAGGCTAAGCGTCGCAA", target_gene="Non-Targeting"),
    ]


@pytest.fixture
def single_guide() -> GuideInput:
    return GuideInput(guide_id="sg_test", spacer_sequence="GAATGGGAACGTACGGGAAA", target_gene="TEST")
