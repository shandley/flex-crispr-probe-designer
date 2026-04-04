"""Guide library validation."""

from flex_crispr_probe_designer.design.dna import gc_content, has_polyt_run, is_valid_dna, max_homopolymer_run
from flex_crispr_probe_designer.models import GuideInput, ValidationResult

EXPECTED_SPACER_LEN = 20
GC_LOW = 0.30
GC_HIGH = 0.70
MAX_HOMOPOLYMER = 5


def validate_guide(guide: GuideInput) -> ValidationResult:
    """Validate a single guide's spacer sequence."""
    errors: list[str] = []
    warnings: list[str] = []

    seq = guide.spacer_sequence

    # Hard errors
    if len(seq) != EXPECTED_SPACER_LEN:
        errors.append(f"Spacer length is {len(seq)}, expected {EXPECTED_SPACER_LEN}")

    if not is_valid_dna(seq):
        errors.append("Spacer contains non-ACGT characters")

    # Warnings
    if is_valid_dna(seq):
        gc = gc_content(seq)
        if gc < GC_LOW or gc > GC_HIGH:
            warnings.append(f"GC content {gc:.1%} outside recommended range ({GC_LOW:.0%}–{GC_HIGH:.0%})")

        if has_polyt_run(seq, min_run=MAX_HOMOPOLYMER):
            warnings.append(f"Contains poly-T run >= {MAX_HOMOPOLYMER} bases (synthesis concern)")

        base, run_len = max_homopolymer_run(seq)
        if run_len >= 6:
            warnings.append(f"Contains {base}x{run_len} homopolymer run")

    return ValidationResult(
        guide_id=guide.guide_id,
        valid=len(errors) == 0,
        warnings=warnings,
        errors=errors,
    )


def validate_guide_library(guides: list[GuideInput]) -> list[ValidationResult]:
    """Validate an entire guide library, including batch-level checks."""
    results = [validate_guide(g) for g in guides]

    # Check for duplicate IDs
    seen_ids: dict[str, int] = {}
    for i, g in enumerate(guides):
        if g.guide_id in seen_ids:
            results[i].errors.append(f"Duplicate guide_id (first seen at index {seen_ids[g.guide_id]})")
            results[i].valid = False
        else:
            seen_ids[g.guide_id] = i

    # Check for duplicate spacers
    seen_spacers: dict[str, str] = {}
    for i, g in enumerate(guides):
        if g.spacer_sequence in seen_spacers:
            results[i].warnings.append(
                f"Duplicate spacer sequence (same as {seen_spacers[g.spacer_sequence]})"
            )
        else:
            seen_spacers[g.spacer_sequence] = g.guide_id

    return results
