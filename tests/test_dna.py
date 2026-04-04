"""Tests for DNA utility functions."""

from __future__ import annotations

from flex_crispr_probe_designer.design.dna import (
    gc_content,
    has_polyt_run,
    is_valid_dna,
    max_homopolymer_run,
    reverse_complement,
)


class TestReverseComplement:
    def test_basic(self) -> None:
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self) -> None:
        assert reverse_complement("AATT") == "AATT"

    def test_single_base(self) -> None:
        assert reverse_complement("A") == "T"

    def test_lowercase_input(self) -> None:
        assert reverse_complement("atcg") == "CGAT"

    def test_known_spacer(self) -> None:
        # GAATGGGAACGTACGGGAAA -> RC = TTTCCCGTACGTTCCCATTC
        assert reverse_complement("GAATGGGAACGTACGGGAAA") == "TTTCCCGTACGTTCCCATTC"


class TestIsValidDna:
    def test_valid(self) -> None:
        assert is_valid_dna("ACGTACGT") is True

    def test_invalid_base(self) -> None:
        assert is_valid_dna("ACGN") is False

    def test_empty(self) -> None:
        assert is_valid_dna("") is True

    def test_lowercase_valid(self) -> None:
        assert is_valid_dna("acgt") is True


class TestGcContent:
    def test_all_gc(self) -> None:
        assert gc_content("GGCC") == 1.0

    def test_all_at(self) -> None:
        assert gc_content("AATT") == 0.0

    def test_half(self) -> None:
        assert gc_content("ACGT") == 0.5

    def test_empty(self) -> None:
        assert gc_content("") == 0.0


class TestMaxHomopolymerRun:
    def test_no_run(self) -> None:
        assert max_homopolymer_run("ACGT") == ("A", 1)

    def test_poly_a(self) -> None:
        assert max_homopolymer_run("AAACGT") == ("A", 3)

    def test_poly_t_at_end(self) -> None:
        assert max_homopolymer_run("ACGTTTTT") == ("T", 5)

    def test_empty(self) -> None:
        assert max_homopolymer_run("") == ("", 0)


class TestHasPolytRun:
    def test_has_run(self) -> None:
        assert has_polyt_run("ACGTTTTTCG") is True

    def test_no_run(self) -> None:
        assert has_polyt_run("ACGTTTTCG") is False

    def test_custom_threshold(self) -> None:
        assert has_polyt_run("ACGTTTTCG", min_run=4) is True
