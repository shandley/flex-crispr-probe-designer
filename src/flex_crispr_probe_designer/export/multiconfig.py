"""Generate Cell Ranger multi config CSVs for Flex v2 CRISPR experiments."""

from itertools import groupby
from pathlib import Path

import pandas as pd

from flex_crispr_probe_designer.models import (
    VALID_COLS,
    VALID_ROWS,
    GemWellConfig,
    MultiConfigParams,
    SampleEntry,
)


def expand_well_range(well_spec: str) -> list[str]:
    """Expand a well specification into a list of well IDs.

    Supports:
      - Single well: "A01"
      - Range: "A01:A12" (same row), "A01:H01" (same column), "A01:B06" (row-major)
      - Comma-separated: "A01,A05,B03"
      - Mixed: "A01:A04,B01,C01:C03"
    """
    wells: list[str] = []
    for part in well_spec.split(","):
        part = part.strip()
        if ":" in part:
            start, end = part.split(":", 1)
            wells.extend(_expand_range(start.strip(), end.strip()))
        else:
            _validate_well(part)
            wells.append(part.upper())
    return wells


def _expand_range(start: str, end: str) -> list[str]:
    """Expand a well range like A01:A12 into individual wells (row-major order)."""
    _validate_well(start)
    _validate_well(end)
    start = start.upper()
    end = end.upper()

    start_row, start_col = start[0], int(start[1:])
    end_row, end_col = end[0], int(end[1:])

    result: list[str] = []
    for row in VALID_ROWS:
        if row < start_row:
            continue
        if row > end_row:
            break
        for col in VALID_COLS:
            if row == start_row and col < start_col:
                continue
            if row == end_row and col > end_col:
                break
            result.append(f"{row}{col:02d}")
    return result


def _validate_well(well: str) -> None:
    """Validate a single well ID like A01."""
    well = well.upper().strip()
    if len(well) < 2 or len(well) > 3:
        raise ValueError(f"Invalid well ID: '{well}' (expected format like A01)")
    row = well[0]
    if row not in VALID_ROWS:
        raise ValueError(f"Invalid well row '{row}' in '{well}' (valid: {VALID_ROWS})")
    try:
        col = int(well[1:])
    except ValueError:
        raise ValueError(f"Invalid well column in '{well}'") from None
    if col not in VALID_COLS:
        raise ValueError(f"Well column {col} out of range 1-12 in '{well}'")


def parse_sample_sheet(path: Path) -> list[SampleEntry]:
    """Parse a sample sheet CSV into SampleEntry objects."""
    df = pd.read_csv(path)
    entries: list[SampleEntry] = []
    for _, row in df.iterrows():
        wells = expand_well_range(str(row["wells"]))
        entries.append(SampleEntry(
            sample_id=str(row["sample_id"]),
            plate=str(row["plate"]),
            wells=wells,
            gem_well=int(row["gem_well"]),
            description=str(row.get("description", "")),
        ))
    return entries


def build_multi_config(samples: list[SampleEntry], params: MultiConfigParams) -> str:
    """Generate a Cell Ranger multi config CSV for one GEM well.

    The config has four sections: [gene-expression], [libraries], [feature], [samples].
    """
    lines: list[str] = []

    # [gene-expression]
    lines.append("[gene-expression]")
    lines.append(f"probe-set,{params.probe_set_path}")
    lines.append(f"create-bam,{'true' if params.create_bam else 'false'}")
    lines.append("")

    # [libraries]
    lines.append("[libraries]")
    lines.append("fastq_id,fastqs,feature_types")
    lines.append(f"{params.gex_fastq_id},{params.gex_fastqs_path},Gene Expression")
    lines.append(f"{params.crispr_fastq_id},{params.crispr_fastqs_path},CRISPR Guide Capture")
    lines.append("")

    # [feature]
    lines.append("[feature]")
    lines.append(f"reference,{params.feature_ref_path}")
    lines.append("")

    # [samples]
    lines.append("[samples]")
    lines.append("sample_id,probe_barcode_ids")
    for s in samples:
        barcode_ids = "|".join(f"{s.plate}-{w}" for w in s.wells)
        lines.append(f"{s.sample_id},{barcode_ids}")

    return "\n".join(lines) + "\n"


def generate_all_configs(
    samples: list[SampleEntry],
    params: MultiConfigParams,
    output_dir: Path,
) -> list[GemWellConfig]:
    """Generate one multi config CSV per GEM well."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Group samples by gem_well
    sorted_samples = sorted(samples, key=lambda s: s.gem_well)
    configs: list[GemWellConfig] = []

    for gem_well, group in groupby(sorted_samples, key=lambda s: s.gem_well):
        well_samples = list(group)
        config_text = build_multi_config(well_samples, params)
        out_path = output_dir / f"multi_config_well{gem_well}.csv"
        out_path.write_text(config_text)
        configs.append(GemWellConfig(
            gem_well=gem_well,
            samples=well_samples,
            config_text=config_text,
            output_path=out_path,
        ))

    return configs
