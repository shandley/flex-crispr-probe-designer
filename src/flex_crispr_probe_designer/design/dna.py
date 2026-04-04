"""DNA sequence utility functions."""

from __future__ import annotations

_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(_COMPLEMENT[b] for b in reversed(seq.upper()))


def is_valid_dna(seq: str) -> bool:
    """Check that a sequence contains only A, C, G, T."""
    return all(b in "ACGT" for b in seq.upper())


def gc_content(seq: str) -> float:
    """Return GC fraction (0.0–1.0) of a DNA sequence."""
    if not seq:
        return 0.0
    upper = seq.upper()
    return sum(1 for b in upper if b in "GC") / len(upper)


def max_homopolymer_run(seq: str) -> tuple[str, int]:
    """Return the base and length of the longest homopolymer run.

    Returns ("", 0) for empty sequences.
    """
    if not seq:
        return ("", 0)
    upper = seq.upper()
    max_base = upper[0]
    max_len = 1
    cur_base = upper[0]
    cur_len = 1
    for b in upper[1:]:
        if b == cur_base:
            cur_len += 1
        else:
            if cur_len > max_len:
                max_base = cur_base
                max_len = cur_len
            cur_base = b
            cur_len = 1
    if cur_len > max_len:
        max_base = cur_base
        max_len = cur_len
    return (max_base, max_len)


def has_polyt_run(seq: str, min_run: int = 5) -> bool:
    """Check if the sequence contains a poly-T run of at least min_run bases."""
    return "T" * min_run in seq.upper()
