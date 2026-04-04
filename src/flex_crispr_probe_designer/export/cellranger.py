"""Generate Cell Ranger feature reference CSV for CRISPR Guide Capture."""

from pathlib import Path

import pandas as pd

from flex_crispr_probe_designer.constants import (
    AT_SCAFFOLD_RC_DNA,
    PREFERRED_SCAFFOLD_RC,
    RHS_LINKER,
)
from flex_crispr_probe_designer.design.dna import reverse_complement
from flex_crispr_probe_designer.models import CustomScaffold, FeatureRefRow, GuideInput, ScaffoldType


def _build_pattern(scaffold: ScaffoldType, custom: CustomScaffold | None = None) -> str:
    """Build the Cell Ranger pattern string for CRISPR Guide Capture.

    The pattern encodes the sequence Cell Ranger expects upstream of the
    spacer barcode in Read 2. For Flex CRISPR, this is:
      scaffold_RC + TAAAC_linker + (BC)

    Verified against 10X's own feature reference from the 1.2M K562 Flex
    CRISPR dataset, where the pattern is:
      GCTATGCTGTTTCCAGCTTAGCTCTTAAAC(BC)
    """
    if scaffold == ScaffoldType.CUSTOM:
        if custom is None:
            raise ValueError("custom scaffold required for custom scaffold type")
        scaffold_rc = reverse_complement(custom.scaffold_sequence)
    elif scaffold == ScaffoldType.AT_SUBSTITUTION:
        scaffold_rc = AT_SCAFFOLD_RC_DNA
    else:
        scaffold_rc = PREFERRED_SCAFFOLD_RC

    return f"{scaffold_rc}{RHS_LINKER}(BC)"


def build_feature_reference(
    guides: list[GuideInput],
    scaffold: ScaffoldType = ScaffoldType.PREFERRED,
    custom: CustomScaffold | None = None,
) -> list[FeatureRefRow]:
    """Build Cell Ranger feature reference rows from a guide library.

    The sequence field contains the protospacer (spacer) sequence.
    The pattern field encodes the scaffold_RC + linker prefix that
    Cell Ranger uses to locate the spacer in Read 2.
    """
    pattern = _build_pattern(scaffold, custom)

    rows: list[FeatureRefRow] = []
    for g in guides:
        target_id = g.target_gene if g.target_gene != "Non-Targeting" else "Non-Targeting"
        target_name = g.target_gene if g.target_gene != "Non-Targeting" else "Non-Targeting"
        rows.append(FeatureRefRow(
            id=g.guide_id,
            name=g.guide_id,
            read="R2",
            pattern=pattern,
            sequence=g.spacer_sequence,
            feature_type="CRISPR Guide Capture",
            target_gene_id=target_id,
            target_gene_name=target_name,
        ))
    return rows


def write_feature_reference(
    guides: list[GuideInput],
    output: Path,
    scaffold: ScaffoldType = ScaffoldType.PREFERRED,
    custom: CustomScaffold | None = None,
) -> Path:
    """Write a Cell Ranger-compatible feature_reference.csv."""
    rows = build_feature_reference(guides, scaffold, custom)
    records = [r.model_dump() for r in rows]
    df = pd.DataFrame(records)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    return output
