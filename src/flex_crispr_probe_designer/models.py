"""Pydantic models for probe design, validation, and experiment planning."""

import enum
from pathlib import Path

from pydantic import BaseModel, Field, field_validator


class ScaffoldType(enum.StrEnum):
    """Built-in guide scaffold variant."""

    PREFERRED = "preferred"
    AT_SUBSTITUTION = "at_substitution"
    CUSTOM = "custom"


class CustomScaffold(BaseModel):
    """User-provided custom guide scaffold sequence.

    The scaffold_sequence should be the actual scaffold DNA (5'→3') as it
    appears in the sgRNA. The tool will compute the reverse complement
    for the LHS probe automatically.

    If idt_notation is provided, it is used as-is for the ordering sheet
    (e.g., with +N LNA modifications). Otherwise, the plain RC is used.
    """

    scaffold_sequence: str = Field(..., description="Scaffold DNA sequence (5'→3') from the sgRNA vector")
    idt_notation: str | None = Field(
        None, description="Optional IDT notation with modifications (e.g., LNA +N). If None, plain RC is used."
    )

    @field_validator("scaffold_sequence")
    @classmethod
    def normalize_scaffold(cls, v: str) -> str:
        v = v.strip().upper()
        valid = set("ACGT")
        if not all(c in valid for c in v):
            raise ValueError(f"Scaffold contains non-ACGT characters: {v}")
        return v


class OPoolScale(enum.StrEnum):
    """IDT oPool synthesis scale."""

    PMOL_50 = "50pmol"
    PMOL_10 = "10pmol"


class GuideInput(BaseModel):
    """A single guide from the input library."""

    guide_id: str
    spacer_sequence: str
    target_gene: str

    @field_validator("spacer_sequence")
    @classmethod
    def normalize_spacer(cls, v: str) -> str:
        return v.strip().upper()


class ValidationResult(BaseModel):
    """Result of validating a single guide."""

    guide_id: str
    valid: bool
    warnings: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)


class LHSProbe(BaseModel):
    """LHS probe — one per experiment (hybridizes to guide scaffold)."""

    scaffold_type: ScaffoldType
    scaffold_label: str = ""  # Human-readable label (e.g., "preferred" or "custom (28 bp)")
    sequence: str  # Pure DNA sequence (no modifications)
    idt_notation: str  # With LNA +N annotations for ordering
    truseq_r2: str
    scaffold_rc: str
    length: int
    has_lna: bool


class RHSProbe(BaseModel):
    """RHS probe — one per guide (recognizes unique spacer)."""

    guide_id: str
    target_gene: str
    spacer_sequence: str
    spacer_rc: str
    sequence: str  # Pure DNA (no modification prefix)
    idt_notation: str  # With /5Phos/ for ordering
    length: int


class ProbeSet(BaseModel):
    """Complete probe set for an experiment."""

    scaffold_type: ScaffoldType
    lhs: LHSProbe
    rhs_probes: list[RHSProbe]
    n_guides: int


class SpikeInPool(BaseModel):
    """Spike-in pool preparation volumes and concentrations."""

    opool_scale: OPoolScale
    scaffold_type: ScaffoldType
    n_hybridizations: int
    idte_vol_ul: float
    rhs_vol_ul: float
    lhs_vol_ul: float
    total_vol_ul: float
    lhs_final_conc_nm: float
    rhs_final_conc_nm: float


class ExperimentPlan(BaseModel):
    """Experiment planning calculations."""

    n_guides: int
    cells_per_guide: int
    min_cells: int
    wta_hybs: int
    bc_oligo_hybs: int
    gem_wells: int
    cells_per_gem_well: int
    cells_hybd_wta: int
    cells_into_fixation: int
    gex_reads_total: int
    crispr_reads_total: int
    total_reads: int


class FeatureRefRow(BaseModel):
    """Single row of Cell Ranger feature reference CSV."""

    id: str
    name: str
    read: str = "R2"
    pattern: str
    sequence: str
    feature_type: str = "CRISPR Guide Capture"
    target_gene_id: str
    target_gene_name: str


class DesignOutput(BaseModel):
    """Paths to all generated output files."""

    ordering_csv: Path
    opool_csv: Path
    feature_reference_csv: Path
    spikein_csv: Path
    experiment_plan_txt: Path


# --- Cell Ranger multi config models ---

VALID_PLATES = {"A", "B", "C", "D"}
VALID_ROWS = "ABCDEFGH"
VALID_COLS = range(1, 13)  # 1-12 for 96-well


class SampleEntry(BaseModel):
    """A sample in the experiment sample sheet."""

    sample_id: str
    plate: str
    wells: list[str]
    gem_well: int
    description: str = ""

    @field_validator("plate")
    @classmethod
    def validate_plate(cls, v: str) -> str:
        v = v.strip().upper()
        if v not in VALID_PLATES:
            raise ValueError(f"Plate must be one of {sorted(VALID_PLATES)}, got '{v}'")
        return v


class MultiConfigParams(BaseModel):
    """Parameters for generating Cell Ranger multi config CSVs."""

    probe_set_path: str
    gex_fastq_id: str = "flex_gex"
    gex_fastqs_path: str
    crispr_fastq_id: str = "flex_cr"
    crispr_fastqs_path: str
    feature_ref_path: str
    create_bam: bool = False


class GemWellConfig(BaseModel):
    """Generated config for a single GEM well."""

    gem_well: int
    samples: list[SampleEntry]
    config_text: str
    output_path: Path | None = None
