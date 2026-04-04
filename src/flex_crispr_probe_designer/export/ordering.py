"""Generate IDT-compatible ordering sheets."""

from pathlib import Path

import pandas as pd

from flex_crispr_probe_designer.models import ProbeSet


def probes_to_ordering_csv(probe_set: ProbeSet, output: Path) -> Path:
    """Write a combined ordering CSV with LHS and all RHS probes.

    Format is compatible with IDT custom oligo ordering.
    """
    rows: list[dict[str, str]] = []

    # LHS probe (single row)
    lhs = probe_set.lhs
    purification = "HPLC" if lhs.has_lna else "Standard Desalted"
    notes = "LNA-modified bases marked with +" if lhs.has_lna else ""
    rows.append({
        "Name": f"LHS_{lhs.scaffold_type.value}",
        "Sequence": lhs.idt_notation,
        "Scale": "25 nmol",
        "Purification": purification,
        "Notes": notes,
    })

    # RHS probes (one per guide)
    for rhs in probe_set.rhs_probes:
        rows.append({
            "Name": f"RHS_{rhs.guide_id}",
            "Sequence": rhs.idt_notation,
            "Scale": "25 nmol",
            "Purification": "Standard Desalted",
            "Notes": f"target={rhs.target_gene}",
        })

    df = pd.DataFrame(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    return output


def probes_to_opool_csv(probe_set: ProbeSet, output: Path) -> Path:
    """Write RHS probes in IDT oPool format.

    oPool format: Name, Sequence columns (all probes in one pool).
    The 50 pmol scale is recommended for most CRISPR workflows.
    """
    rows: list[dict[str, str]] = []
    pool_name = f"CRISPR_RHS_{probe_set.scaffold_type.value}"

    for rhs in probe_set.rhs_probes:
        rows.append({
            "Pool Name": pool_name,
            "Sequence Name": rhs.guide_id,
            "Sequence": rhs.idt_notation,
        })

    df = pd.DataFrame(rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    return output
