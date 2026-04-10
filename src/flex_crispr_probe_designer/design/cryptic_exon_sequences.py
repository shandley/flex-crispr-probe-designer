"""Pre-fetched hg38 genomic sequences for cryptic exon junction probe design.

Sequences fetched from UCSC Genome Browser API (hg38).
Each entry provides 30 bp from each side of the splice junctions.

For probe design, we use the inner 25 bp from each side to create
50 bp junction probes (matching Flex v2 WTA custom probe format per CG000839).
"""

from flex_crispr_probe_designer.design.cryptic_exons import (
    CrypticExonProbe,
    design_cryptic_exon_probes_from_sequences,
)

# === STMN2 (chr8, + strand) ===
# Coordinates from Baughn et al. 2023, verified against UCSC hg38

# Junction 1: Exon 1 end (79,611,185-79,611,214) | CE start (79,616,822-79,616,851)
STMN2_UPSTREAM_EXON_3PRIME = "ACATCCCTACAATGGCTAAAACAGCAATGG"  # 30 bp
STMN2_CE_5PRIME = "GACTCGGCAGAAGACCTTCGAGAGAAAGGT"  # 30 bp

# Junction 2: CE end (79,617,019-79,617,048) | Exon 2 start (79,622,174-79,622,203)
STMN2_CE_3PRIME = "GAAGAAATTAAAACAAAAGATTGCTGTCTC"  # 30 bp
STMN2_DOWNSTREAM_EXON_5PRIME = "GAATTCTTAAATTTATTAATCATTTCTATA"  # 30 bp

# === UNC13A (chr19, - strand) ===
# Sequences are reverse-complemented to mRNA sense orientation
# Coordinates from Ma et al. 2022 Nature, Brown et al. 2022 Nature

# Junction 1: Upstream exon 3' | CE 5' (mRNA sense, RC of genomic)
UNC13A_UPSTREAM_EXON_3PRIME = "AGCCATGACGTGAGACTGCAGGCTGGGGTG"  # 30 bp
UNC13A_CE_5PRIME = "GAGTAGATAAAAGGATGGATGGAGAGATGG"  # 30 bp

# Junction 2: CE 3' | Downstream exon 5' (mRNA sense)
UNC13A_CE_3PRIME = "CTGCCTGGGTTTCCTGGAAAGAACTCTTAT"  # 30 bp
UNC13A_DOWNSTREAM_EXON_5PRIME = "TGGAAGCATCATTTGTACCTGGGAGGTTGA"  # 30 bp


# === HDGFL2 (chr19, + strand) ===
# Coordinates from Seddighi et al. 2024 Sci Transl Med
# NOTE: CE boundaries are approximate -- the 1,867 bp intron region needs
# RNA-seq confirmation. Confidence: MEDIUM.

# Junction 1: Exon 6 end (4,491,806-4,491,835) | CE start (4,491,836-4,491,865)
HDGFL2_UPSTREAM_EXON_3PRIME = "AGGGGCCCTCTGGGGGGACGGAAAAAAAAG"  # 30 bp
HDGFL2_CE_5PRIME = "GTAGCGTGCACTTGACTTTGTTTCCCATGC"  # 30 bp

# Junction 2: CE end (4,493,673-4,493,702) | Exon 7 start (4,493,703-4,493,732)
HDGFL2_CE_3PRIME = "ATGCTCACGCCTGTCCCTGCTTCTCCCCAG"  # 30 bp
HDGFL2_DOWNSTREAM_EXON_5PRIME = "AAGGCGCCATCAGCCTCCGACTCCGACTCC"  # 30 bp

# === KALRN (chr3, + strand) ===
# Coordinates from Humphrey et al. 2023 Acta Neuropathol

# Junction 1: Upstream exon end (124,701,225-124,701,254) | CE start (124,701,255-124,701,284)
KALRN_UPSTREAM_EXON_3PRIME = "TCATTGCCAAGGAATTGACACATCACACAA"  # 30 bp
KALRN_CE_5PRIME = "GGTAATACTGCCAAAGAATCTACTCAGGCA"  # 30 bp

# Junction 2: CE end (124,702,009-124,702,038) | Downstream exon start (124,702,039-124,702,068)
KALRN_CE_3PRIME = "ATTGTAAATGGATTCTCTTTCTTCCGAAGA"  # 30 bp
KALRN_DOWNSTREAM_EXON_5PRIME = "TGCTGCTGCTGATGGTGCCACCATTTCTTG"  # 30 bp

# === SYT7 (chr11, - strand) ===
# Sequences are reverse-complemented to mRNA sense orientation
# Coordinates from Humphrey et al. 2023 Acta Neuropathol

# Junction 1: Upstream exon 3' | CE 5' (mRNA sense, RC of genomic)
SYT7_UPSTREAM_EXON_3PRIME = "TGATCTAGACAGAGACTTTTGGAATAACAA"  # 30 bp
SYT7_CE_5PRIME = "GTTTTCTCCCCTCTCTCTCTGCCCCTCCAG"  # 30 bp

# Junction 2: CE 3' | Downstream exon 5' (mRNA sense)
SYT7_CE_3PRIME = "GTATGTTTGCACATGTGTGCGTGTGCGTGT"  # 30 bp
SYT7_DOWNSTREAM_EXON_5PRIME = "CGGCCCTCGAGTGCGTGCGTGTGTGCATTG"  # 30 bp


def generate_stmn2_probes() -> list[CrypticExonProbe]:
    """Generate junction probes for STMN2 cryptic exon.

    Returns 3 probes:
    - inclusion_5prime: Exon1-CE junction (detects CE inclusion)
    - inclusion_3prime: CE-Exon2 junction (detects CE inclusion)
    - skipping: Exon1-Exon2 junction (detects normal splice, no CE)
    """
    return design_cryptic_exon_probes_from_sequences(
        gene="STMN2",
        upstream_exon_3prime=STMN2_UPSTREAM_EXON_3PRIME,
        ce_5prime=STMN2_CE_5PRIME,
        ce_3prime=STMN2_CE_3PRIME,
        downstream_exon_5prime=STMN2_DOWNSTREAM_EXON_5PRIME,
    )


def generate_unc13a_probes() -> list[CrypticExonProbe]:
    """Generate junction probes for UNC13A cryptic exon.

    Returns 3 probes (in mRNA sense orientation, already RC'd for - strand).
    """
    return design_cryptic_exon_probes_from_sequences(
        gene="UNC13A",
        upstream_exon_3prime=UNC13A_UPSTREAM_EXON_3PRIME,
        ce_5prime=UNC13A_CE_5PRIME,
        ce_3prime=UNC13A_CE_3PRIME,
        downstream_exon_5prime=UNC13A_DOWNSTREAM_EXON_5PRIME,
    )


def generate_hdgfl2_probes() -> list[CrypticExonProbe]:
    """Generate junction probes for HDGFL2 cryptic exon.

    NOTE: CE boundaries are approximate (1,867 bp intron, MEDIUM confidence).
    Exact boundaries need RNA-seq confirmation per Seddighi et al. 2024.

    Returns 3 probes: 5' inclusion, 3' inclusion, and skipping control.
    """
    return design_cryptic_exon_probes_from_sequences(
        gene="HDGFL2",
        upstream_exon_3prime=HDGFL2_UPSTREAM_EXON_3PRIME,
        ce_5prime=HDGFL2_CE_5PRIME,
        ce_3prime=HDGFL2_CE_3PRIME,
        downstream_exon_5prime=HDGFL2_DOWNSTREAM_EXON_5PRIME,
    )


def generate_kalrn_probes() -> list[CrypticExonProbe]:
    """Generate junction probes for KALRN cryptic exon.

    Returns 3 probes: 5' inclusion, 3' inclusion, and skipping control.
    """
    return design_cryptic_exon_probes_from_sequences(
        gene="KALRN",
        upstream_exon_3prime=KALRN_UPSTREAM_EXON_3PRIME,
        ce_5prime=KALRN_CE_5PRIME,
        ce_3prime=KALRN_CE_3PRIME,
        downstream_exon_5prime=KALRN_DOWNSTREAM_EXON_5PRIME,
    )


def generate_syt7_probes() -> list[CrypticExonProbe]:
    """Generate junction probes for SYT7 cryptic exon.

    Returns 3 probes (in mRNA sense orientation, already RC'd for - strand).
    """
    return design_cryptic_exon_probes_from_sequences(
        gene="SYT7",
        upstream_exon_3prime=SYT7_UPSTREAM_EXON_3PRIME,
        ce_5prime=SYT7_CE_5PRIME,
        ce_3prime=SYT7_CE_3PRIME,
        downstream_exon_5prime=SYT7_DOWNSTREAM_EXON_5PRIME,
    )


def generate_all_cryptic_exon_probes() -> dict[str, list[CrypticExonProbe]]:
    """Generate probes for all cryptic exons with available sequences."""
    return {
        "STMN2": generate_stmn2_probes(),
        "UNC13A": generate_unc13a_probes(),
        "HDGFL2": generate_hdgfl2_probes(),
        "KALRN": generate_kalrn_probes(),
        "SYT7": generate_syt7_probes(),
    }
