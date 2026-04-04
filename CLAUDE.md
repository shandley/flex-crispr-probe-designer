# flex-crispr-probe-designer

CRISPR guide capture probe designer for 10X Genomics Flex v2.

## Quick reference

```bash
uv run pytest tests/              # Run tests (83 tests)
uv run ruff check src/ tests/     # Lint
uv run pyright src/               # Type check (strict mode)
uv run flex-crispr --help         # CLI usage
```

## CLI commands

```bash
# Design probes from a guide library
flex-crispr design guides.csv -o output/ --scaffold preferred
flex-crispr design guides.csv -o output/ --custom-scaffold "AGAGCTATGCTGGAAACAGCATAGC"

# Plan experiment (cells, reagents, sequencing)
flex-crispr plan 5000

# Validate a guide library
flex-crispr validate guides.csv

# Generate Cell Ranger multi config CSVs
flex-crispr config sample_sheet.csv --probe-set /path/to/probes.csv \
  --gex-fastqs /data/gex --crispr-fastqs /data/crispr \
  --feature-ref /data/feature_reference.csv -o output/
```

## Architecture

- `constants.py` — All biological constants from 10X SAM001062 technical note
- `models.py` — Pydantic models (strict validation throughout)
- `design/dna.py` — DNA utilities (RC, GC, validation)
- `design/probes.py` — LHS/RHS probe generation (preferred, A>T, custom scaffolds)
- `validate/guides.py` — Guide library validation
- `export/ordering.py` — IDT ordering sheets (individual + oPool)
- `export/cellranger.py` — Cell Ranger feature reference CSV
- `export/multiconfig.py` — Cell Ranger multi config CSV generation
- `export/spikein.py` — Spike-in pool preparation tables
- `plan/experiment.py` — Experiment planner (cells, hybs, GEM wells, seq depth)
- `cli.py` — Typer CLI with rich output

## Key technical notes

- Probe sequences derived from SAM001062 Rev A (CRISPR with Flex v2)
- Feature reference uses `5P(BC)` pattern — verify against Cell Ranger v10.0 docs before production use
- Constants should be verified against CG000814 technical note before ordering oligos
- Supports human and mouse only (10X probe set constraint)
