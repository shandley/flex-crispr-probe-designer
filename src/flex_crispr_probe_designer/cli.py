"""CLI for flex-crispr probe designer."""

from pathlib import Path

import pandas as pd
import typer
from rich.console import Console
from rich.table import Table

from flex_crispr_probe_designer.design.cryptic_exons import get_cryptic_exon, list_cryptic_exons
from flex_crispr_probe_designer.design.probes import build_probe_set
from flex_crispr_probe_designer.export.cellranger import write_feature_reference
from flex_crispr_probe_designer.export.multiconfig import generate_all_configs, parse_sample_sheet
from flex_crispr_probe_designer.export.ordering import probes_to_opool_csv, probes_to_ordering_csv
from flex_crispr_probe_designer.export.spikein import calculate_spikein_pool, write_spikein_csv
from flex_crispr_probe_designer.models import (
    CustomScaffold,
    ExperimentPlan,
    GuideInput,
    MultiConfigParams,
    OPoolScale,
    ScaffoldType,
)
from flex_crispr_probe_designer.plan.experiment import plan_experiment
from flex_crispr_probe_designer.validate.guides import validate_guide_library

app = typer.Typer(name="flex-crispr", help="CRISPR guide capture probe designer for 10X Genomics Flex v2")
console = Console()

# Column name aliases for flexible CSV parsing
_ID_ALIASES = {"guide_id", "id", "name", "sgrna_id", "sgrna", "guide"}
_SEQ_ALIASES = {"spacer_sequence", "spacer", "sequence", "protospacer", "seq"}
_GENE_ALIASES = {"target_gene", "gene", "target", "target_gene_name", "gene_name"}


def _find_column(columns: list[str], aliases: set[str]) -> str | None:
    """Find a column name matching any alias (case-insensitive)."""
    lower_map = {c.lower().strip(): c for c in columns}
    for alias in aliases:
        if alias.lower() in lower_map:
            return lower_map[alias.lower()]
    return None


def _parse_guide_library(input_path: Path) -> list[GuideInput]:
    """Parse a CSV/TSV guide library with flexible column naming."""
    sep = "\t" if input_path.suffix in (".tsv", ".txt") else ","
    df = pd.read_csv(input_path, sep=sep)
    cols: list[str] = [str(c) for c in df.columns]

    id_col = _find_column(cols, _ID_ALIASES)
    seq_col = _find_column(cols, _SEQ_ALIASES)
    gene_col = _find_column(cols, _GENE_ALIASES)

    if not id_col or not seq_col or not gene_col:
        missing: list[str] = []
        if not id_col:
            missing.append(f"guide ID (tried: {', '.join(sorted(_ID_ALIASES))})")
        if not seq_col:
            missing.append(f"spacer sequence (tried: {', '.join(sorted(_SEQ_ALIASES))})")
        if not gene_col:
            missing.append(f"target gene (tried: {', '.join(sorted(_GENE_ALIASES))})")
        console.print(f"[red]Could not find required columns: {'; '.join(missing)}[/red]")
        console.print(f"[dim]Found columns: {', '.join(cols)}[/dim]")
        raise typer.Exit(1)

    guides: list[GuideInput] = []
    for _, row in df.iterrows():
        guides.append(GuideInput(
            guide_id=str(row[id_col]),
            spacer_sequence=str(row[seq_col]),
            target_gene=str(row[gene_col]),
        ))
    return guides


@app.command()
def design(
    input_path: Path = typer.Argument(..., help="Guide library CSV/TSV"),
    output_dir: Path = typer.Option("output", "-o", "--output", help="Output directory"),
    scaffold: ScaffoldType = typer.Option(ScaffoldType.PREFERRED, "--scaffold", "-s", help="Scaffold type"),
    custom_scaffold: str | None = typer.Option(
        None, "--custom-scaffold",
        help="Custom scaffold DNA sequence (5'->3' as in sgRNA). Implies --scaffold custom.",
    ),
    custom_scaffold_idt: str | None = typer.Option(
        None, "--custom-scaffold-idt",
        help="IDT notation for custom scaffold RC (e.g., with +N LNA mods). Used for ordering sheet.",
    ),
    opool_scale: OPoolScale = typer.Option(OPoolScale.PMOL_50, "--opool-scale", help="oPool synthesis scale"),
) -> None:
    """Design probes from a guide library and generate all output files."""
    # Handle custom scaffold
    custom: CustomScaffold | None = None
    if custom_scaffold:
        scaffold = ScaffoldType.CUSTOM
        custom = CustomScaffold(scaffold_sequence=custom_scaffold, idt_notation=custom_scaffold_idt)

    console.print(f"[bold]Reading guide library:[/bold] {input_path}")
    guides = _parse_guide_library(input_path)
    console.print(f"  Loaded {len(guides)} guides")

    # Validate
    results = validate_guide_library(guides)
    n_errors = sum(1 for r in results if not r.valid)
    n_warnings = sum(1 for r in results if r.warnings)

    if n_errors > 0:
        console.print(f"[red]  {n_errors} guides with errors:[/red]")
        for r in results:
            if not r.valid:
                for e in r.errors:
                    console.print(f"    [red]{r.guide_id}: {e}[/red]")
        raise typer.Exit(1)

    if n_warnings > 0:
        console.print(f"[yellow]  {n_warnings} guides with warnings[/yellow]")

    # Build probes
    scaffold_label = custom.scaffold_sequence if custom else scaffold.value
    console.print(f"[bold]Building probes:[/bold] scaffold={scaffold_label}")
    probe_set = build_probe_set(guides, scaffold, custom)

    # Export all files
    output_dir.mkdir(parents=True, exist_ok=True)

    ordering_path = probes_to_ordering_csv(probe_set, output_dir / "probe_ordering.csv")
    opool_path = probes_to_opool_csv(probe_set, output_dir / "rhs_opool.csv")
    feat_ref_path = write_feature_reference(guides, output_dir / "feature_reference.csv", scaffold, custom)

    pool = calculate_spikein_pool(scaffold, opool_scale)
    spikein_path = write_spikein_csv(pool, output_dir / "spikein_pool.csv")

    # Experiment plan
    exp_plan = plan_experiment(len(guides))
    plan_path = output_dir / "experiment_plan.txt"
    plan_path.write_text(_format_plan(exp_plan, scaffold))

    # Summary
    console.print("\n[bold green]Output files:[/bold green]")
    for p in [ordering_path, opool_path, feat_ref_path, spikein_path, plan_path]:
        console.print(f"  {p}")

    console.print("\n[bold]Probe summary:[/bold]")
    console.print(f"  LHS probe: {probe_set.lhs.length} bp ({probe_set.lhs.scaffold_label} scaffold)")
    console.print(f"  RHS probes: {probe_set.n_guides} x {probe_set.rhs_probes[0].length} bp")
    if probe_set.lhs.has_lna:
        console.print("  [yellow]Note: LHS contains LNA modifications (requires HPLC purification)[/yellow]")
    if scaffold == ScaffoldType.CUSTOM:
        target_len = len(probe_set.lhs.scaffold_rc) + 20  # scaffold RC + spacer portion
        if target_len != 50:
            console.print(
                f"  [yellow]Note: targeting portion is {target_len} bp (not 50). "
                f"Read 2 must be >= {target_len} cycles.[/yellow]"
            )


@app.command()
def plan(
    n_guides: int = typer.Argument(..., help="Number of guides in the library"),
    cells_per_guide: int = typer.Option(200, "--cells-per-guide", "-c", help="Target cells per guide"),
) -> None:
    """Calculate experiment requirements for a given number of guides."""
    exp = plan_experiment(n_guides, cells_per_guide)

    table = Table(title=f"Experiment Plan: {n_guides:,} guides ({cells_per_guide} cells/guide)")
    table.add_column("Parameter", style="bold")
    table.add_column("Value", justify="right")

    table.add_row("Min cells needed", f"{exp.min_cells:,}")
    table.add_row("WTA hybridizations", str(exp.wta_hybs))
    table.add_row("Barcode Oligo hybridizations", str(exp.bc_oligo_hybs))
    table.add_row("GEM wells", str(exp.gem_wells))
    table.add_row("Cells per GEM well", f"{exp.cells_per_gem_well:,}")
    table.add_row("Cells into WTA hyb", f"{exp.cells_hybd_wta:,}")
    table.add_row("Cells into fixation", f"{exp.cells_into_fixation:,}")
    table.add_row("", "")
    table.add_row("GEX reads (10k/cell)", f"{exp.gex_reads_total:,}")
    table.add_row("CRISPR reads (2k/cell)", f"{exp.crispr_reads_total:,}")
    table.add_row("Total reads", f"{exp.total_reads:,}")

    console.print(table)


@app.command()
def validate(
    input_path: Path = typer.Argument(..., help="Guide library CSV/TSV"),
) -> None:
    """Validate a guide library without generating probes."""
    guides = _parse_guide_library(input_path)
    results = validate_guide_library(guides)

    n_valid = sum(1 for r in results if r.valid)
    n_errors = sum(1 for r in results if not r.valid)
    n_warnings = sum(1 for r in results if r.warnings)

    console.print(f"[bold]Validated {len(guides)} guides[/bold]")
    console.print(
        f"  [green]{n_valid} valid[/green]  [red]{n_errors} errors[/red]  [yellow]{n_warnings} warnings[/yellow]"
    )

    for r in results:
        if r.errors or r.warnings:
            console.print(f"\n  [bold]{r.guide_id}:[/bold]")
            for e in r.errors:
                console.print(f"    [red]ERROR: {e}[/red]")
            for w in r.warnings:
                console.print(f"    [yellow]WARN: {w}[/yellow]")


@app.command()
def config(
    sample_sheet: Path = typer.Argument(..., help="Sample sheet CSV (sample_id, plate, wells, gem_well)"),
    output_dir: Path = typer.Option("output", "-o", "--output", help="Output directory"),
    probe_set: str = typer.Option(..., "--probe-set", help="Path to probe set CSV"),
    gex_fastqs: str = typer.Option(..., "--gex-fastqs", help="Path to Gene Expression FASTQ directory"),
    crispr_fastqs: str = typer.Option(..., "--crispr-fastqs", help="Path to CRISPR Guide Capture FASTQ directory"),
    feature_ref: str = typer.Option(..., "--feature-ref", help="Path to feature reference CSV"),
    gex_fastq_id: str = typer.Option("flex_gex", "--gex-fastq-id", help="FASTQ ID for GEX libraries"),
    crispr_fastq_id: str = typer.Option("flex_cr", "--crispr-fastq-id", help="FASTQ ID for CRISPR libraries"),
    create_bam: bool = typer.Option(False, "--create-bam", help="Generate BAM files (default: false)"),
) -> None:
    """Generate Cell Ranger multi config CSVs from a sample sheet."""
    console.print(f"[bold]Reading sample sheet:[/bold] {sample_sheet}")
    samples = parse_sample_sheet(sample_sheet)

    n_wells = len({s.gem_well for s in samples})
    total_barcodes = sum(len(s.wells) for s in samples)
    console.print(f"  {len(samples)} samples, {total_barcodes} barcodes, {n_wells} GEM wells")

    params = MultiConfigParams(
        probe_set_path=probe_set,
        gex_fastq_id=gex_fastq_id,
        gex_fastqs_path=gex_fastqs,
        crispr_fastq_id=crispr_fastq_id,
        crispr_fastqs_path=crispr_fastqs,
        feature_ref_path=feature_ref,
        create_bam=create_bam,
    )

    configs = generate_all_configs(samples, params, output_dir)

    console.print(f"\n[bold green]Generated {len(configs)} config files:[/bold green]")
    for c in configs:
        n_samples = len(c.samples)
        n_barcodes = sum(len(s.wells) for s in c.samples)
        console.print(f"  GEM well {c.gem_well}: {n_samples} samples, {n_barcodes} barcodes → {c.output_path}")


@app.command(name="cryptic-exons")
def cryptic_exons(
    gene: str | None = typer.Argument(None, help="Gene name to look up (e.g., STMN2). Omit to list all."),
) -> None:
    """List known TDP-43 cryptic exons and their coordinates."""
    if gene is None:
        # List all
        table = Table(title="Known TDP-43 Cryptic Exons (hg38)")
        table.add_column("Gene", style="bold")
        table.add_column("Chr")
        table.add_column("Start")
        table.add_column("End")
        table.add_column("Length")
        table.add_column("Strand")
        table.add_column("Location")
        table.add_column("Confidence")

        for ce in list_cryptic_exons():
            table.add_row(
                ce.gene, ce.chrom,
                f"{ce.ce_start:,}", f"{ce.ce_end:,}",
                f"{ce.ce_length:,} bp", ce.strand,
                ce.intron_location, ce.confidence,
            )
        console.print(table)
        console.print(
            "\n[dim]Use [bold]flex-crispr cryptic-exons GENE[/bold] for details. "
            "Probe design requires genomic sequences from these coordinates.[/dim]"
        )
    else:
        ce = get_cryptic_exon(gene)
        if ce is None:
            console.print(f"[red]Gene '{gene}' not found in cryptic exon catalog.[/red]")
            console.print("[dim]Available: " + ", ".join(c.gene for c in list_cryptic_exons()) + "[/dim]")
            raise typer.Exit(1)

        console.print(f"\n[bold]{ce.gene}[/bold] ({ce.gene_id})")
        console.print(f"  Chromosome:  {ce.chrom}:{ce.ce_start:,}-{ce.ce_end:,} ({ce.strand})")
        console.print(f"  CE length:   {ce.ce_length:,} bp")
        console.print(f"  Location:    {ce.intron_location}")
        console.print(f"  Mechanism:   {ce.mechanism}")
        console.print(f"  Confidence:  {ce.confidence}")
        console.print(f"  Reference:   {ce.key_reference}")
        if ce.notes:
            console.print(f"  Notes:       {ce.notes}")

        console.print("\n[bold]Probe design:[/bold]")
        console.print("  To design junction probes, extract 25 bp sequences from each side")
        console.print("  of the splice junctions using the coordinates above, then use:")
        console.print(f"    design_cryptic_exon_probes_from_sequences('{ce.gene}', ...)")
        console.print(f"\n  Upstream exon ends at:     {ce.chrom}:{ce.upstream_exon_end:,}")
        console.print(f"  Cryptic exon:              {ce.chrom}:{ce.ce_start:,}-{ce.ce_end:,}")
        console.print(f"  Downstream exon starts at: {ce.chrom}:{ce.downstream_exon_start:,}")


def _format_plan(exp: ExperimentPlan, scaffold: ScaffoldType) -> str:
    """Format experiment plan as plain text."""
    from flex_crispr_probe_designer.constants import KIT_PN_16_SAMPLE, SEQ_READ1_CONFIG

    lines = [
        "Flex v2 CRISPR Experiment Plan",
        "=" * 40,
        f"Guides:                  {exp.n_guides:,}",
        f"Cells per guide:         {exp.cells_per_guide}",
        f"Scaffold:                {scaffold.value}",
        "",
        "Cell Requirements",
        "-" * 40,
        f"Min cells needed:        {exp.min_cells:,}",
        f"Cells into WTA hyb:      {exp.cells_hybd_wta:,}",
        f"Cells into fixation:     {exp.cells_into_fixation:,}",
        "",
        "Reagents",
        "-" * 40,
        f"WTA hybridizations:      {exp.wta_hybs}",
        f"Barcode Oligo hybs:      {exp.bc_oligo_hybs}",
        f"GEM wells:               {exp.gem_wells}",
        f"Cells per GEM well:      {exp.cells_per_gem_well:,}",
        "",
        "Sequencing (Read 1 config)",
        "-" * 40,
        f"Read 1:                  {SEQ_READ1_CONFIG['read1']} cycles",
        f"i7 Index:                {SEQ_READ1_CONFIG['i7_index']} cycles",
        f"i5 Index:                {SEQ_READ1_CONFIG['i5_index']} cycles",
        f"Read 2:                  {SEQ_READ1_CONFIG['read2']} cycles",
        f"GEX reads (10k/cell):    {exp.gex_reads_total:,}",
        f"CRISPR reads (2k/cell):  {exp.crispr_reads_total:,}",
        f"Total reads:             {exp.total_reads:,}",
        "",
        f"Suggested kit: GEM-X Flex v2 Human, 16 samples ({KIT_PN_16_SAMPLE})",
    ]
    return "\n".join(lines) + "\n"
