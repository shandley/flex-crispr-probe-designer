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


def generate_all_cryptic_exon_probes() -> dict[str, list[CrypticExonProbe]]:
    """Generate probes for all cryptic exons with available sequences."""
    return {
        "STMN2": generate_stmn2_probes(),
        "UNC13A": generate_unc13a_probes(),
    }
