"""Biological and protocol constants from 10X Genomics SAM001062 (Flex v2 CRISPR)."""

from __future__ import annotations

# === Probe sequence components ===

# TruSeq Read 2 handle (34 bp) — 5' end of LHS probe
TRUSEQ_R2_HANDLE = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

# Guide scaffold reverse complements for LHS probe
# Preferred scaffold RC (25 bp)
PREFERRED_SCAFFOLD_RC = "GCTATGCTGTTTCCAGCTTAGCTCT"
# A>T substitution scaffold RC with LNA positions marked by + prefix (28 bp)
# The + before a base denotes a Locked Nucleic Acid modification
AT_SCAFFOLD_RC_IDT = "+C+TTGCTATGC+TGT+TT+CC+AGCATAGCTCT"
# Same sequence without LNA annotations (for pure DNA operations)
AT_SCAFFOLD_RC_DNA = "CTTGCTATGCTGTTTCCAGCATAGCTCT"

# RHS probe components
RHS_FIVE_PHOS = "/5Phos/"  # IDT 5' phosphorylation notation
RHS_LINKER = "TAAAC"  # 5 bp linker at 5' of RHS after phosphate
RHS_PCONSTANT_SEQ = "CCCATATAAGAAA"  # 13 bp partial Constant Sequence at 3' of RHS

# === Spike-in pool volumes (µL) for 16 hybridizations ===

SPIKEIN_VOLUMES: dict[str, dict[str, float]] = {
    "50pmol": {"idte_ph8": 39.6, "rhs_stock": 2.2, "lhs_8um": 2.2, "total": 44.0},
    "10pmol": {"idte_ph8": 30.8, "rhs_stock": 11.0, "lhs_8um": 2.2, "total": 44.0},
}

# Final spike-in concentrations (nM)
LHS_CONC_PREFERRED_NM = 400
LHS_CONC_AT_SUBST_NM = 1200
RHS_CONC_PER_PROBE_NM = 40

# To achieve 1200 nM LHS for A>T scaffold: increase LHS volume to 6.6 µL,
# reduce IDTE proportionally
LHS_VOL_AT_SCAFFOLD_UL = 6.6

# === Experiment planning table (SAM001062 slide 25) ===
# Assumes 1 biological sample subpooled for WTA hyb setup, 200 cells/guide

EXPERIMENT_PLANNING: list[dict[str, int]] = [
    {"guides": 100, "min_cells": 20_000, "wta_hybs": 1, "bc_oligo_hybs": 1, "gem_wells": 1,
     "cells_per_gem_well": 20_000, "cells_hybd_wta": 42_000, "cells_into_fixation": 60_000},
    {"guides": 400, "min_cells": 80_000, "wta_hybs": 1, "bc_oligo_hybs": 4, "gem_wells": 1,
     "cells_per_gem_well": 80_000, "cells_hybd_wta": 165_000, "cells_into_fixation": 240_000},
    {"guides": 1_200, "min_cells": 240_000, "wta_hybs": 1, "bc_oligo_hybs": 12, "gem_wells": 1,
     "cells_per_gem_well": 240_000, "cells_hybd_wta": 500_000, "cells_into_fixation": 710_000},
    {"guides": 2_400, "min_cells": 480_000, "wta_hybs": 2, "bc_oligo_hybs": 24, "gem_wells": 1,
     "cells_per_gem_well": 480_000, "cells_hybd_wta": 1_000_000, "cells_into_fixation": 1_400_000},
    {"guides": 4_800, "min_cells": 960_000, "wta_hybs": 4, "bc_oligo_hybs": 48, "gem_wells": 1,
     "cells_per_gem_well": 960_000, "cells_hybd_wta": 2_000_000, "cells_into_fixation": 2_800_000},
    {"guides": 9_600, "min_cells": 1_920_000, "wta_hybs": 8, "bc_oligo_hybs": 48, "gem_wells": 2,
     "cells_per_gem_well": 960_000, "cells_hybd_wta": 4_000_000, "cells_into_fixation": 5_700_000},
    {"guides": 14_000, "min_cells": 2_800_000, "wta_hybs": 12, "bc_oligo_hybs": 48, "gem_wells": 3,
     "cells_per_gem_well": 933_333, "cells_hybd_wta": 5_800_000, "cells_into_fixation": 8_300_000},
    {"guides": 19_000, "min_cells": 3_800_000, "wta_hybs": 16, "bc_oligo_hybs": 48, "gem_wells": 4,
     "cells_per_gem_well": 950_000, "cells_hybd_wta": 7_900_000, "cells_into_fixation": 11_300_000},
]

# === Sequencing parameters ===

# Read 1 Sequencing Configuration (preferred for Flex v2 CRISPR)
SEQ_READ1_CONFIG = {"read1": 54, "i7_index": 10, "i5_index": 10, "read2": 50}

# Minimum read pairs per cell
GEX_MIN_READS_PER_CELL = 10_000
CRISPR_REC_READS_PER_CELL = 2_000
PROTEIN_MIN_READS_PER_CELL = 5_000

# Default cells per targeting sgRNA
DEFAULT_CELLS_PER_GUIDE = 200

# === User-supplied TruSeq primers (for Feature PCR) ===

TRUSEQ_R1_FWD_PRIMER = "CTACACGACGCTCTTCCGATCT"
TRUSEQ_R2_REV_PRIMER = "GTGACTGGAGTTCAGACGTGTG"

# === Kit part numbers ===

KIT_PN_16_SAMPLE = "PN-1000927"  # GEM-X Flex v2 Human, 16 samples
KIT_PN_96_SAMPLE = "PN-1000928"  # GEM-X Flex v2 Human, 96 samples
KIT_PN_DUAL_INDEX_TT = "PN-1000215"  # Dual Index Kit TT Set A, 96 rxns
KIT_PN_CORE_REAGENT = "PN-1000937"  # GEM-X Flex Core Reagent Bundle v2, 96 rxns
