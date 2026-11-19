"""
Command-line interface.
=======================

A single Typer application exposes one subcommand per pipeline stage. Snakemake
rules call these subcommands with explicit paths, so the CLI never has to read
the dataset declarations embedded in ``workflow/datasets/``.

Input
-----
- Paths supplied by the caller and ``config/config.yaml``.

Output
------
- Exit code ``0`` on success, ``1`` on a caught fatal error.
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
import hashlib
from pathlib import Path

import click
from loguru import logger
import pandas as pd
import typer

from src.afdb import (
    CIF_ENTRY_COLUMNS,
    PAE_INDEX_COLUMNS,
    join_entry_table,
    read_cif_entries,
    read_pae_index,
)
from src.dna import (
    build_genome_bed,
    build_nucleotide_bed,
    map_domains_to_nucleotides,
    organize_domains,
)
from src.external import build_ecod_table, build_ted_table
from src.files import (
    AppConfig,
    load_config,
    load_gene_table,
    read_tsv,
    require_columns,
    setup_logger,
    write_tsv,
)
from src.pae import build_domain_table, run_pae_to_domains, run_pae_to_domains_from_registry
from src.plddt import build_plddt_table, build_plddt_table_from_registry
from src.registry import RegistryPaths, build_registry
from src.segments import (
    annotate_table,
    build_confident_table,
    build_reformated_table,
    build_reformated_table_from_registry,
)
from src.unify import (
    combine_segments,
    unify_segments as build_unified_segments,
    unify_segments_from_registry,
)
from src.uniprot import (
    RESOLUTION_INPUT_COLUMNS,
    TRIANGULATION_INPUT_COLUMNS,
    UniprotCol,
    build_resolution,
    build_triangulation,
    fetch_entries,
)

# =============================================================================
# 2. CONSTANTS
# =============================================================================
DEFAULT_CONFIG = Path("config/config.yaml")
SHA256_SIZE_LIMIT = 50 * 1024 * 1024
app = typer.Typer(add_completion=False, no_args_is_help=True)


# =============================================================================
# 3. IMPERATIVE SHELL (CLI Commands)
# =============================================================================
def _sha256(path: Path) -> str:
    """Hash small resources; large dumps are identified by size only."""
    if path.stat().st_size > SHA256_SIZE_LIMIT:
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resource_manifest(config: AppConfig) -> pd.DataFrame:
    """Describe every configured resource with size and checksum."""
    entries = {
        "pombase_gene_table": config.resources.pombase_gene_table,
        "pombase_peptide_stats": config.resources.pombase_peptide_stats,
        "pombase_gff3": config.resources.pombase_gff3,
        "pombase_features_gene_table": config.resources.pombase_features_gene_table,
        "pombase_features_peptide_stats": config.resources.pombase_features_peptide_stats,
        "ted_dump": config.resources.ted_dump,
        "ted_segmentation_dump": config.resources.ted_segmentation_dump,
        "ecod_dump": config.resources.ecod_dump,
    }
    records = [
        {
            "name": name,
            "path": str(path),
            "exists": path.exists(),
            "size_bytes": path.stat().st_size if path.exists() else 0,
            "sha256": _sha256(path) if path.exists() else "",
        }
        for name, path in entries.items()
    ]
    return pd.DataFrame.from_records(records)


@app.command("fetch-resources")
def fetch_resources(
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
    output: Path = typer.Option(Path("resources/manifest.tsv"), help="Manifest output path."),
) -> None:
    """Verify every configured resource and write a provenance manifest."""
    manifest = _resource_manifest(load_config(config))
    write_tsv(manifest, output)
    logger.info(f"Resource manifest written to {output}")


@app.command("afdb-cif-entries")
def afdb_cif_entries(
    afdb_dir: Path = typer.Option(..., help="AFdb v6 CIF directory."),
    output: Path = typer.Option(..., help="Output CIF entry table."),
) -> None:
    """Read every AFdb CIF header into the entry inventory."""
    entries = read_cif_entries(afdb_dir)
    write_tsv(entries, output)
    logger.info(f"AFdb CIF entries written to {output} ({len(entries):,} rows)")


@app.command("afdb-pae-index")
def afdb_pae_index(
    pae_dir: Path = typer.Option(..., help="AFdb PAE JSON directory."),
    output: Path = typer.Option(..., help="Output PAE index table."),
) -> None:
    """Index the AFdb PAE matrices by accession."""
    index = read_pae_index(pae_dir)
    write_tsv(index, output)
    logger.info(f"AFdb PAE index written to {output} ({len(index):,} rows)")


@app.command("afdb-entry-table")
def afdb_entry_table(
    cif_entries: Path = typer.Option(..., help="AFdb CIF entry table."),
    pae_index: Path = typer.Option(..., help="AFdb PAE index table."),
    output: Path = typer.Option(..., help="Output AFdb entry table."),
) -> None:
    """Join PAE availability onto the AFdb CIF inventory."""
    entries = read_tsv(cif_entries, dtype=str)
    index = read_tsv(pae_index, dtype=str)
    require_columns(entries, CIF_ENTRY_COLUMNS, cif_entries)
    require_columns(index, PAE_INDEX_COLUMNS, pae_index)

    table = join_entry_table(entries, index)
    write_tsv(table, output)
    logger.info(f"AFdb entry table written to {output} ({len(table):,} rows)")


@app.command("refresh-uniprot-triangulation")
def refresh_uniprot_triangulation(
    nonidentical_pairs: Path = typer.Option(..., help="Non-identical pair table."),
    output: Path = typer.Option(..., help="Refresh candidate output."),
) -> None:
    """Rebuild the UniProt triangulation snapshot from the current pairs."""
    pairs = read_tsv(nonidentical_pairs, dtype=str)
    require_columns(pairs, TRIANGULATION_INPUT_COLUMNS, nonidentical_pairs)

    accessions = pairs[UniprotCol.UNIPROT_ID].tolist()
    entries = fetch_entries(accessions)

    table = build_triangulation(pairs, entries)
    write_tsv(table, output)
    logger.info(f"UniProt triangulation refresh written to {output} ({len(table):,} rows)")


@app.command("refresh-uniprot-resolution")
def refresh_uniprot_resolution(
    afdb_only: Path = typer.Option(..., help="AFdb-only entry table."),
    gene_table: Path = typer.Option(..., help="PomBase gene table."),
    output: Path = typer.Option(..., help="Refresh candidate output."),
) -> None:
    """Rebuild the UniProt accession-resolution snapshot for AFdb-only entries."""
    queries = read_tsv(afdb_only, dtype=str)
    require_columns(queries, RESOLUTION_INPUT_COLUMNS, afdb_only)
    is_superseded = queries["match_status"] == "superseded_by_sequence_match"
    queries = queries.loc[~is_superseded].reset_index(drop=True)
    pombase_accessions = set(load_gene_table(gene_table)[UniprotCol.UNIPROT_ID].dropna())

    accessions = queries[UniprotCol.UNIPROT_ID].tolist()
    entries = fetch_entries(accessions)

    table = build_resolution(queries, entries, pombase_accessions)
    write_tsv(table, output)
    logger.info(f"UniProt resolution refresh written to {output} ({len(table):,} rows)")


@app.command("build-registry")
def build_registry_command(
    merged_table: Path = typer.Option(..., help="Comparison merged table."),
    output: Path = typer.Option(..., help="Registry output table."),
    summary_output: Path = typer.Option(..., help="Selection summary output table."),
    stats_output: Path = typer.Option(..., help="Downstream statistics output table."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Build the per-gene structure registry."""
    paths = RegistryPaths(
        merged_table=merged_table,
        output=output,
        summary_output=summary_output,
        stats_output=stats_output,
    )
    build_registry(paths, load_config(config))


@app.command("registry-pae-to-domains")
def registry_pae_to_domains(
    registry: Path = typer.Option(..., help="Gene structure registry."),
    output_dir: Path = typer.Option(..., help="Directory for .domain files."),
    jobs: int = typer.Option(1, min=1, help="Parallel worker processes."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Cluster the PAE matrix selected for every registry gene."""
    run_pae_to_domains_from_registry(registry, output_dir, jobs, load_config(config))


@app.command("registry-domains-to-table")
def registry_domains_to_table(
    domains_dir: Path = typer.Option(..., help="Directory of .domain files."),
    output: Path = typer.Option(..., help="Output PAE domain table."),
) -> None:
    """Collect registry .domain files into the PAE domain table."""
    build_domain_table(domains_dir, output, id_column="gene_systematic_id")


@app.command("registry-extract-plddt")
def registry_extract_plddt(
    registry: Path = typer.Option(..., help="Gene structure registry."),
    output: Path = typer.Option(..., help="Output pLDDT table."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Extract pLDDT for the structure selected for every registry gene."""
    build_plddt_table_from_registry(registry, output, load_config(config))


@app.command("registry-merge-domains")
def registry_merge_domains(
    pae_table: Path = typer.Option(..., help="PAE domain table."),
    plddt_table: Path = typer.Option(..., help="pLDDT table."),
    registry: Path = typer.Option(..., help="Gene structure registry."),
    output: Path = typer.Option(..., help="Output reformatted domain table."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Merge registry PAE clusters with pLDDT confidence."""
    build_reformated_table_from_registry(
        pae_table, plddt_table, registry, load_config(config).parameters, output
    )


@app.command("registry-unify-segments")
def registry_unify_segments(
    pae_table: Path = typer.Option(..., help="Annotated PAE table."),
    ted_table: Path = typer.Option(..., help="PomBase TED table."),
    ecod_table: Path = typer.Option(..., help="PomBase ECOD table."),
    registry: Path = typer.Option(..., help="Gene structure registry."),
    output: Path = typer.Option(..., help="Output unified segment table."),
    dataset: str = typer.Option(..., help="Dataset identifier."),
    pae_version: str = typer.Option(..., help="PAE model version label."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Build the unified segment table from registry-selected domains."""
    app_config = load_config(config)
    unify_segments_from_registry(
        pae_table,
        ted_table,
        ecod_table,
        registry,
        output,
        dataset=dataset,
        pae_source_version=pae_version,
        ted_version=app_config.resources.ted_release,
        ecod_version=app_config.resources.ecod_release,
    )


@app.command("pae-to-domains")
def pae_to_domains(
    pae_dir: Path = typer.Option(..., help="Directory of PAE JSON files."),
    output_dir: Path = typer.Option(..., help="Directory for .domain files."),
    jobs: int = typer.Option(1, min=1, help="Parallel worker processes."),
    protein_subset: Path | None = typer.Option(None, help="Optional one-ID-per-line subset."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Cluster PAE matrices into domain files."""
    run_pae_to_domains(pae_dir, output_dir, load_config(config).parameters, jobs, protein_subset)


@app.command("domains-to-table")
def domains_to_table(
    domains_dir: Path = typer.Option(..., help="Directory of .domain files."),
    output: Path = typer.Option(..., help="Output PAE domain table."),
) -> None:
    """Collect .domain files into the PAE domain table."""
    build_domain_table(domains_dir, output)


@app.command("extract-plddt")
def extract_plddt(
    structure_dir: Path = typer.Option(..., help="Directory of PDB/mmCIF structures."),
    output: Path = typer.Option(..., help="Output pLDDT table."),
    protein_subset: Path | None = typer.Option(None, help="Optional one-ID-per-line subset."),
) -> None:
    """Extract per-residue pLDDT values from structures."""
    build_plddt_table(structure_dir, output, protein_subset)


@app.command("merge-domains")
def merge_domains(
    pae_table: Path = typer.Option(..., help="PAE domain table."),
    plddt_table: Path = typer.Option(..., help="pLDDT table."),
    gene_table: Path = typer.Option(..., help="PomBase gene table."),
    output: Path = typer.Option(..., help="Output reformatted domain table."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Merge PAE clusters with pLDDT confidence."""
    build_reformated_table(
        pae_table, plddt_table, gene_table, load_config(config).parameters, output
    )


@app.command("assign-type")
def assign_type(
    input_table: Path = typer.Option(..., "--input", help="Reformatted domain table."),
    output: Path = typer.Option(..., help="Output annotated domain table."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Assign domain location and domain/disordered type."""
    annotate_table(input_table, load_config(config).parameters, output)


@app.command("confident-domains")
def confident_domains(
    plddt_table: Path = typer.Option(..., help="pLDDT table."),
    gene_table: Path = typer.Option(..., help="PomBase gene table."),
    output: Path = typer.Option(..., help="Output confident-domain table."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Run the independent 3D-AF-Surfer style confident-domain caller."""
    build_confident_table(plddt_table, gene_table, load_config(config).parameters, output)


@app.command("dna-bed")
def dna_bed(
    gff: Path = typer.Option(..., help="Genome GFF3."),
    peptide_stats: Path = typer.Option(..., help="PomBase PeptideStats table."),
    output: Path = typer.Option(..., help="Output per-nucleotide table."),
) -> None:
    """Create the per-nucleotide coding-region table."""
    build_nucleotide_bed(gff, peptide_stats, output)


@app.command("dna-map")
def dna_map(
    nucleotide_bed: Path = typer.Option(..., help="Per-nucleotide table."),
    annotated: Path = typer.Option(..., help="Annotated protein domain table."),
    output: Path = typer.Option(..., help="Output nucleotide-to-domain table."),
) -> None:
    """Map protein domains onto nucleotides."""
    map_domains_to_nucleotides(nucleotide_bed, annotated, output)


@app.command("dna-organize")
def dna_organize(
    input_table: Path = typer.Option(..., "--input", help="Nucleotide-to-domain table."),
    output: Path = typer.Option(..., help="Output merged DNA domain table."),
) -> None:
    """Merge nucleotide-level assignments into genomic intervals."""
    organize_domains(input_table, output)


@app.command("dna-genome-bed")
def dna_genome_bed(
    input_table: Path = typer.Option(..., "--input", help="Merged DNA domain table."),
    output: Path = typer.Option(..., help="Output BED9 table."),
) -> None:
    """Expand merged genomic intervals into BED9."""
    build_genome_bed(input_table, output)


@app.command("ted-pombe")
def ted_pombe(
    ted_dump: Path = typer.Option(..., help="TED detailed dump."),
    gene_table: Path = typer.Option(..., help="PomBase gene table."),
    output: Path = typer.Option(..., help="Output PomBase TED table."),
) -> None:
    """Extract the PomBase subset of the TED dump."""
    build_ted_table(ted_dump, gene_table, output)


@app.command("ecod-pombe")
def ecod_pombe(
    ecod_dump: Path = typer.Option(..., help="ECOD domain dump."),
    gene_table: Path = typer.Option(..., help="PomBase gene table."),
    output: Path = typer.Option(..., help="Output PomBase ECOD table."),
) -> None:
    """Extract the PomBase subset of the ECOD dump."""
    build_ecod_table(ecod_dump, gene_table, output)


@app.command("unify-segments")
def unify_segments(
    pae_table: Path = typer.Option(..., help="Annotated PAE table."),
    ted_table: Path = typer.Option(..., help="PomBase TED table."),
    ecod_table: Path = typer.Option(..., help="PomBase ECOD table."),
    gene_table: Path = typer.Option(..., help="PomBase gene table for annotation."),
    stats_table: Path = typer.Option(..., help="PomBase PeptideStats table."),
    output: Path = typer.Option(..., help="Output unified segment table."),
    dataset: str = typer.Option(..., help="Dataset identifier."),
    pae_version: str = typer.Option(..., help="PAE model version label."),
    config: Path = typer.Option(DEFAULT_CONFIG, help="Project configuration file."),
) -> None:
    """Build the unified segment table for one dataset."""
    app_config = load_config(config)
    build_unified_segments(
        pae_table,
        ted_table,
        ecod_table,
        gene_table,
        stats_table,
        output,
        dataset=dataset,
        pae_source_version=pae_version,
        ted_version=app_config.resources.ted_release,
        ecod_version=app_config.resources.ecod_release,
    )


@app.command("combine-segments")
def combine_segments_command(
    inputs: list[Path] = typer.Option(..., "--input", help="Per-dataset segment tables."),
    output: Path = typer.Option(..., help="Combined output table."),
) -> None:
    """Concatenate per-dataset segment tables."""
    combine_segments(inputs, output)


# =============================================================================
# 4. ENTRY POINT
# =============================================================================
def main() -> int:
    """Dispatch the Typer application and propagate the exit code."""
    setup_logger()
    try:
        app(prog_name="afdb-domains", standalone_mode=False)
    except click.exceptions.Exit as exit_signal:
        return int(exit_signal.exit_code)
    except click.ClickException as error:
        error.show()
        return error.exit_code
    except Exception as error:
        logger.exception(f"Fatal execution error: {error}")
        return 1
    return 0
