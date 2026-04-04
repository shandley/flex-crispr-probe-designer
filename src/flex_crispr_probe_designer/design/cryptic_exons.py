"""Custom WTA probe design for known cryptic exons.

Designs probe pairs for Flex v2 spike-in where one probe sits in the
cryptic exon and the partner sits in the flanking constitutive exon.
Enables detection of cryptic exon inclusion at single-cell resolution.

Coordinates are GRCh38/hg38. These cryptic exons are NOT in standard
GENCODE annotations — they are only included upon TDP-43 loss of function.
"""

from pydantic import BaseModel

from flex_crispr_probe_designer.design.dna import gc_content, reverse_complement


class CrypticExon(BaseModel):
    """A known cryptic exon with coordinates and probe design info."""

    gene: str
    gene_id: str
    chrom: str
    ce_start: int  # 0-based, hg38
    ce_end: int
    strand: str  # + or -
    ce_length: int
    intron_location: str  # e.g., "intron 1 (between exons 1-2)"
    upstream_exon_end: int  # Last base of upstream constitutive exon
    downstream_exon_start: int  # First base of downstream constitutive exon
    mechanism: str
    confidence: str  # HIGH, MEDIUM, LOW
    key_reference: str
    notes: str = ""


class CrypticExonProbe(BaseModel):
    """A designed probe pair targeting a cryptic exon junction."""

    gene: str
    probe_name: str
    ce_sequence: str  # Sequence from cryptic exon (25 bp)
    flanking_sequence: str  # Sequence from constitutive exon (25 bp)
    junction_probe: str  # Combined 50 bp probe spanning the junction
    junction_probe_rc: str
    gc_content: float
    probe_type: str  # "inclusion" (CE present) or "skipping" (CE absent)
    description: str


# === Known cryptic exon catalog ===
# Sources: Humphrey et al. 2023, Seddighi et al. 2024, Baughn et al. 2023,
#          Ma et al. 2022, Brown et al. 2022

CRYPTIC_EXON_CATALOG: list[CrypticExon] = [
    CrypticExon(
        gene="STMN2",
        gene_id="ENSG00000104435",
        chrom="chr8",
        ce_start=79_616_822,
        ce_end=79_617_048,
        strand="+",
        ce_length=227,
        intron_location="intron 1 (between exons 1-2)",
        upstream_exon_end=79_611_214,
        downstream_exon_start=79_622_174,
        mechanism="TDP-43 loss → cryptic exon inclusion → premature stop + polyA",
        confidence="HIGH",
        key_reference="Baughn et al. 2023 Science; Klim et al. 2019 Nat Neurosci",
        notes="Most established ALS/FTD biomarker. ASO therapeutics in clinical trials. "
              "CE contains stop codon + polyA signal 203 bp apart. "
              "TDP-43 binds GU-rich motif with 3x GUGUGU hexamers ~5.6 kb into intron.",
    ),
    CrypticExon(
        gene="UNC13A",
        gene_id="ENSG00000130477",
        chrom="chr19",
        ce_start=17_642_414,
        ce_end=17_642_541,
        strand="-",
        ce_length=128,
        intron_location="intron 20 (between exons 20-21)",
        upstream_exon_end=17_641_401,
        downstream_exon_start=17_644_272,
        mechanism="TDP-43 loss → cryptic exon inclusion → NMD",
        confidence="HIGH",
        key_reference="Ma et al. 2022 Nature; Brown et al. 2022 Nature",
        notes="ALS/FTD GWAS risk SNP rs12973192 is within the CE itself. "
              "128 bp and 178 bp variants share same 3' end. Not conserved in mice.",
    ),
    CrypticExon(
        gene="HDGFL2",
        gene_id="ENSG00000173821",
        chrom="chr19",
        ce_start=4_491_836,
        ce_end=4_493_702,
        strand="+",
        ce_length=1867,
        intron_location="intron 6 (between exons 6-7)",
        upstream_exon_end=4_491_835,
        downstream_exon_start=4_493_703,
        mechanism="TDP-43 loss → in-frame cryptic exon → 46 aa cryptic peptide",
        confidence="MEDIUM",
        key_reference="Seddighi et al. 2024 Sci Transl Med",
        notes="Only known in-frame CE producing stable cryptic peptide detectable in CSF. "
              "Monoclonal antibody exists (CST #95895). "
              "Exact CE boundaries within the 1,867 bp intron need RNA-seq confirmation.",
    ),
    CrypticExon(
        gene="KALRN",
        gene_id="ENSG00000166263",
        chrom="chr3",
        ce_start=124_701_255,
        ce_end=124_702_038,
        strand="+",
        ce_length=784,
        intron_location="near 3' end of gene",
        upstream_exon_end=124_701_254,
        downstream_exon_start=124_702_039,
        mechanism="TDP-43 loss → cryptic exon inclusion",
        confidence="HIGH",
        key_reference="Humphrey et al. 2023 Acta Neuropathol",
        notes="Very large gene (~693 kb). CE 3' end coincides with downstream canonical exon boundary.",
    ),
    CrypticExon(
        gene="SYT7",
        gene_id="ENSG00000011347",
        chrom="chr11",
        ce_start=61_547_309,
        ce_end=61_547_529,
        strand="-",
        ce_length=221,
        intron_location="intron 3-4",
        upstream_exon_end=61_547_308,
        downstream_exon_start=61_547_530,
        mechanism="TDP-43 loss → cryptic exon inclusion",
        confidence="HIGH",
        key_reference="Humphrey et al. 2023 Acta Neuropathol",
        notes="CE immediately adjacent to exon 4 (1 bp gap). Reverse strand.",
    ),
]


def list_cryptic_exons() -> list[CrypticExon]:
    """Return the catalog of known cryptic exons."""
    return CRYPTIC_EXON_CATALOG.copy()


def get_cryptic_exon(gene: str) -> CrypticExon | None:
    """Look up a cryptic exon by gene name (case-insensitive)."""
    gene_upper = gene.upper()
    for ce in CRYPTIC_EXON_CATALOG:
        if ce.gene.upper() == gene_upper:
            return ce
    return None


def design_junction_probe(
    ce_sequence_at_junction: str,
    flanking_sequence_at_junction: str,
    gene: str,
    probe_type: str = "inclusion",
) -> CrypticExonProbe:
    """Design a probe spanning a cryptic exon-constitutive exon junction.

    For an "inclusion" probe:
    - ce_sequence_at_junction: 25 bp from the cryptic exon side of the junction
    - flanking_sequence_at_junction: 25 bp from the constitutive exon side

    The combined 50 bp probe will only hybridize when the cryptic exon is included.

    For Flex v2 WTA custom probes, this 50 bp sequence is split into
    LHS (25 bp) + RHS (25 bp) probe pairs per CG000839.
    """
    junction = flanking_sequence_at_junction + ce_sequence_at_junction
    junction_rc = reverse_complement(junction)
    gc = gc_content(junction)

    if probe_type == "inclusion":
        desc = f"{gene} cryptic exon inclusion junction probe"
    else:
        desc = f"{gene} exon skipping junction probe (CE absent)"

    return CrypticExonProbe(
        gene=gene,
        probe_name=f"{gene}_CE_{probe_type}",
        ce_sequence=ce_sequence_at_junction,
        flanking_sequence=flanking_sequence_at_junction,
        junction_probe=junction,
        junction_probe_rc=junction_rc,
        gc_content=round(gc, 3),
        probe_type=probe_type,
        description=desc,
    )


def design_cryptic_exon_probes_from_sequences(
    gene: str,
    upstream_exon_3prime: str,
    ce_5prime: str,
    ce_3prime: str,
    downstream_exon_5prime: str,
) -> list[CrypticExonProbe]:
    """Design probes for a cryptic exon given flanking sequences.

    Requires 25 bp sequences from each side of two junctions:
    1. Upstream exon → CE (5' junction): upstream_exon_3prime + ce_5prime
    2. CE → downstream exon (3' junction): ce_3prime + downstream_exon_5prime

    Also designs a "skipping" probe for the normal splice (CE absent):
    3. Upstream exon → downstream exon: upstream_exon_3prime + downstream_exon_5prime

    Returns 3 probes: 5' inclusion, 3' inclusion, and skipping control.
    """
    probes: list[CrypticExonProbe] = []

    # 5' inclusion junction: upstream_exon | CE_start
    probes.append(design_junction_probe(
        ce_sequence_at_junction=ce_5prime[:25],
        flanking_sequence_at_junction=upstream_exon_3prime[-25:],
        gene=gene,
        probe_type="inclusion_5prime",
    ))

    # 3' inclusion junction: CE_end | downstream_exon
    probes.append(design_junction_probe(
        ce_sequence_at_junction=ce_3prime[-25:],
        flanking_sequence_at_junction=downstream_exon_5prime[:25],
        gene=gene,
        probe_type="inclusion_3prime",
    ))

    # Skipping junction: upstream_exon | downstream_exon (no CE)
    probes.append(design_junction_probe(
        ce_sequence_at_junction=downstream_exon_5prime[:25],
        flanking_sequence_at_junction=upstream_exon_3prime[-25:],
        gene=gene,
        probe_type="skipping",
    ))

    return probes
