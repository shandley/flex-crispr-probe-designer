"""Tests for guide validation."""

from flex_crispr_probe_designer.models import GuideInput
from flex_crispr_probe_designer.validate.guides import validate_guide, validate_guide_library


class TestValidateGuide:
    def test_valid_guide(self) -> None:
        g = GuideInput(guide_id="sg1", spacer_sequence="GAATGGGAACGTACGGGAAA", target_gene="BRCA1")
        r = validate_guide(g)
        assert r.valid is True
        assert len(r.errors) == 0

    def test_wrong_length(self) -> None:
        g = GuideInput(guide_id="sg1", spacer_sequence="ACGT", target_gene="X")
        r = validate_guide(g)
        assert r.valid is False
        assert any("length" in e.lower() for e in r.errors)

    def test_invalid_bases(self) -> None:
        g = GuideInput(guide_id="sg1", spacer_sequence="GAATGGGAACGTACGGGAAN", target_gene="X")
        r = validate_guide(g)
        assert r.valid is False
        assert any("non-ACGT" in e for e in r.errors)

    def test_low_gc_warning(self) -> None:
        # 20 bp, all A/T = 0% GC
        g = GuideInput(guide_id="sg1", spacer_sequence="AATTAATTAATTAATTAATT", target_gene="X")
        r = validate_guide(g)
        assert r.valid is True  # Warning, not error
        assert any("GC content" in w for w in r.warnings)

    def test_polyt_warning(self) -> None:
        g = GuideInput(guide_id="sg1", spacer_sequence="ACGTTTTTTCGACGACGACG", target_gene="X")
        r = validate_guide(g)
        assert any("poly-T" in w for w in r.warnings)


class TestValidateLibrary:
    def test_duplicate_id(self) -> None:
        guides = [
            GuideInput(guide_id="sg1", spacer_sequence="GAATGGGAACGTACGGGAAA", target_gene="A"),
            GuideInput(guide_id="sg1", spacer_sequence="TCCTCAAGGAGCTTTATTGA", target_gene="B"),
        ]
        results = validate_guide_library(guides)
        assert results[1].valid is False
        assert any("Duplicate guide_id" in e for e in results[1].errors)

    def test_duplicate_spacer_warning(self) -> None:
        guides = [
            GuideInput(guide_id="sg1", spacer_sequence="GAATGGGAACGTACGGGAAA", target_gene="A"),
            GuideInput(guide_id="sg2", spacer_sequence="GAATGGGAACGTACGGGAAA", target_gene="B"),
        ]
        results = validate_guide_library(guides)
        assert results[1].valid is True  # Warning, not error
        assert any("Duplicate spacer" in w for w in results[1].warnings)

    def test_valid_library(self, sample_guides: list[GuideInput]) -> None:
        results = validate_guide_library(sample_guides)
        assert all(r.valid for r in results)
