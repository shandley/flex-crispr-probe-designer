# flex-crispr-probe-designer

CRISPR guide capture probe designer for 10X Genomics Chromium Flex v2.

Takes a guide library CSV and generates everything needed for a Flex v2 CRISPR screen: probe ordering sheets, Cell Ranger configs, spike-in pool tables, and experiment plans.

## Installation

```bash
# Requires Python 3.12+ and uv
git clone https://github.com/shandley/flex-crispr-probe-designer.git
cd flex-crispr-probe-designer
uv sync
```

## Quick Start

```bash
# Design probes from a guide library
flex-crispr design guides.csv -o output/

# Plan an experiment
flex-crispr plan 5000

# Validate a guide library
flex-crispr validate guides.csv

# Generate Cell Ranger multi configs
flex-crispr config sample_sheet.csv \
  --probe-set /path/to/probe_set.csv \
  --gex-fastqs /data/gex \
  --crispr-fastqs /data/crispr \
  --feature-ref output/feature_reference.csv \
  -o output/

# List known TDP-43 cryptic exons
flex-crispr cryptic-exons
flex-crispr cryptic-exons STMN2
```

## Commands

| Command | Description |
|---------|-------------|
| `design` | Generate probe ordering sheets, feature reference, spike-in tables |
| `plan` | Calculate cells, reagents, sequencing for a given library size |
| `validate` | Check a guide library for errors and warnings |
| `config` | Generate Cell Ranger multi config CSVs from a sample sheet |
| `cryptic-exons` | Look up known TDP-43 cryptic exons with hg38 coordinates |

## Input Format

Guide library CSV with flexible column naming:

```csv
guide_id,spacer_sequence,target_gene
sg_BRCA1_1,GAATGGGAACGTACGGGAAA,BRCA1
sg_NT_1,ACGGAGGCTAAGCGTCGCAA,Non-Targeting
```

Accepted column names: `guide_id`/`id`/`name`, `spacer_sequence`/`spacer`/`protospacer`, `target_gene`/`gene`/`target`.

## Output Files

| File | Description |
|------|-------------|
| `probe_ordering.csv` | IDT-ready ordering sheet (LHS + all RHS with `/5Phos/` and LNA annotations) |
| `rhs_opool.csv` | IDT oPool format for bulk RHS probe ordering |
| `feature_reference.csv` | Cell Ranger `cellranger multi` compatible feature reference |
| `spikein_pool.csv` | Reagent volumes for spike-in pool preparation |
| `experiment_plan.txt` | Cells, WTA hybs, barcode oligos, GEM wells, sequencing depth |

## Features

- **Scaffold support**: Preferred, A>T substitution (with LNA), and custom scaffolds
- **Verified patterns**: Cell Ranger feature reference `pattern` field validated against real 10X Flex CRISPR data
- **Flexible parsing**: Auto-detects CSV column names from common aliases
- **Experiment planner**: Interpolates from 10X's recommended planning table (SAM001062)
- **Cryptic exon catalog**: 5 well-characterized TDP-43 cryptic exons with hg38 coordinates and junction probe design

## Custom Scaffold

```bash
# Provide your vector's scaffold sequence (5'->3' as in the sgRNA)
flex-crispr design guides.csv --custom-scaffold "AGAGCTATGCTGGAAACAGCATAGC"

# With LNA modifications for ordering
flex-crispr design guides.csv \
  --custom-scaffold "AGAGCTATGCTGGAAACAGCATAGC" \
  --custom-scaffold-idt "+GCTATGCTGTTTC+CAGCTTAGCTCT"
```

## Development

```bash
uv run pytest tests/      # 107 tests
uv run ruff check src/    # Lint
uv run pyright src/       # Type check (strict)
```

## Reference

Probe sequences and protocol constants derived from 10X Genomics SAM001062 Rev A ("CRISPR with Flex v2: Custom Probe Design and Workflow Overview") and technical note CG000814.

## License

MIT
