"""Core probe generation logic for Flex v2 CRISPR guide capture."""

from flex_crispr_probe_designer.constants import (
    AT_SCAFFOLD_RC_DNA,
    AT_SCAFFOLD_RC_IDT,
    PREFERRED_SCAFFOLD_RC,
    RHS_FIVE_PHOS,
    RHS_LINKER,
    RHS_PCONSTANT_SEQ,
    TRUSEQ_R2_HANDLE,
)
from flex_crispr_probe_designer.design.dna import reverse_complement
from flex_crispr_probe_designer.models import (
    CustomScaffold,
    GuideInput,
    LHSProbe,
    ProbeSet,
    RHSProbe,
    ScaffoldType,
)


def build_lhs_probe(
    scaffold: ScaffoldType,
    custom: CustomScaffold | None = None,
) -> LHSProbe:
    """Build the LHS probe for a given scaffold type.

    One LHS probe is shared across all guides in an experiment.
    It hybridizes to the constant guide scaffold region.

    For custom scaffolds, the scaffold_sequence is reverse-complemented
    to create the LHS binding region. If idt_notation is provided on the
    CustomScaffold, it is used for the ordering sheet.
    """
    if scaffold == ScaffoldType.CUSTOM:
        if custom is None:
            raise ValueError("custom scaffold must be provided when scaffold type is 'custom'")
        scaffold_rc = reverse_complement(custom.scaffold_sequence)
        dna_seq = TRUSEQ_R2_HANDLE + scaffold_rc
        if custom.idt_notation:
            idt_seq = TRUSEQ_R2_HANDLE + custom.idt_notation
            has_lna = "+" in custom.idt_notation
        else:
            idt_seq = dna_seq
            has_lna = False
        label = f"custom ({len(custom.scaffold_sequence)} bp)"
    elif scaffold == ScaffoldType.PREFERRED:
        scaffold_rc = PREFERRED_SCAFFOLD_RC
        dna_seq = TRUSEQ_R2_HANDLE + scaffold_rc
        idt_seq = dna_seq
        has_lna = False
        label = "preferred"
    else:
        scaffold_rc = AT_SCAFFOLD_RC_DNA
        dna_seq = TRUSEQ_R2_HANDLE + scaffold_rc
        idt_seq = TRUSEQ_R2_HANDLE + AT_SCAFFOLD_RC_IDT
        has_lna = True
        label = "at_substitution"

    return LHSProbe(
        scaffold_type=scaffold,
        scaffold_label=label,
        sequence=dna_seq,
        idt_notation=idt_seq,
        truseq_r2=TRUSEQ_R2_HANDLE,
        scaffold_rc=scaffold_rc,
        length=len(dna_seq),
        has_lna=has_lna,
    )


def build_rhs_probe(guide: GuideInput) -> RHSProbe:
    """Build an RHS probe for a single guide.

    Each guide gets a unique RHS probe that recognizes its spacer sequence.
    Structure: /5Phos/TAAAC-[spacer_RC]-CCCATATAAGAAA
    """
    spacer_rc = reverse_complement(guide.spacer_sequence)
    dna_seq = RHS_LINKER + spacer_rc + RHS_PCONSTANT_SEQ
    idt_seq = RHS_FIVE_PHOS + dna_seq

    return RHSProbe(
        guide_id=guide.guide_id,
        target_gene=guide.target_gene,
        spacer_sequence=guide.spacer_sequence,
        spacer_rc=spacer_rc,
        sequence=dna_seq,
        idt_notation=idt_seq,
        length=len(dna_seq),
    )


def build_probe_set(
    guides: list[GuideInput],
    scaffold: ScaffoldType,
    custom: CustomScaffold | None = None,
) -> ProbeSet:
    """Build a complete probe set from a guide library."""
    lhs = build_lhs_probe(scaffold, custom)
    rhs_probes = [build_rhs_probe(g) for g in guides]
    return ProbeSet(
        scaffold_type=scaffold,
        lhs=lhs,
        rhs_probes=rhs_probes,
        n_guides=len(guides),
    )
