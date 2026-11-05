"""
AFdb v6 vs PomBase vs DY47073 Protein Annotation Consistency
============================================================

Compare the proteins modelled in AlphaFold DB release v6 for
Schizosaccharomyces pombe (proteome UP000002485) with the current PomBase
protein-coding gene annotation.

The comparison has three levels:

1. Set membership by UniProt accession: shared, AFdb-only, PomBase-only.
2. Protein length agreement for shared proteins.
3. Sequence agreement, classified as identical, substitution-only,
   indel-only, or both (alignment-derived).

AFdb sequences are read from the modelled CIF entity sequence
(``_entity_poly.pdbx_seq_one_letter_code_can``); PomBase sequences are read
from ``peptide.fa`` and truncated at the first stop codon. AFdb-only entries
are additionally matched to PomBase peptides by exact sequence, which
separates obsolete/secondary UniProt accessions from genuinely AFdb-specific
proteins. The PAE directory is checked for one-to-one accession pairing.

When a DY47073 protein FASTA is supplied, the report adds a third annotation
layer. DY47073 records are resolved to PomBase gene IDs, collapsed per gene,
and classified against both PomBase and the primary AFdb sequence: all three
identical, DY47073 agrees with PomBase only, DY47073 agrees with AFdb only,
DY47073 differs from both, or DY47073 has no record. DY47073-only proteins
(mostly PomBase 2026 pseudogenes) and AFdb-only entries that match a DY47073
protein are reported separately.

Input
-----
- AFdb directory containing ``AF-<accession>-F<n>-model_v6.cif.gz`` files.
- AFdb PAE directory containing ``AF-<accession>-F<n>-predicted_aligned_error_v6.json``.
- PomBase release directory containing ``gene_IDs_names_products.tsv`` and
  ``peptide.fa[.gz]``.
- Optional DY47073 protein FASTA (one protein record per gene or isoform).
- Optional UniProt accession resolution TSV for the AFdb-only entries.
- Optional AFdb/PomBase/current-UniProt triangulation TSV for shared proteins.

Output
------
- ``afdb_vs_pombase.xlsx``: multi-sheet workbook with the merged table and
  per-category sheets.
- ``report.md``: human-readable difference report.
- ``summary_counts.tsv``: headline counts.
- ``merged_table.tsv``: unified shared/AFdb-only/PomBase-only table.
- ``shared_comparison.tsv``: one row per ID-matched shared gene pair.
- ``shared_by_sequence.tsv``: shared pairs matched by exact sequence.
- ``nonidentical_pairs.tsv``: shared pairs with non-identical sequences.
- ``afdb_only.tsv``: AFdb-only accessions with exact-match evidence.
- ``pombase_only.tsv``: PomBase-only genes with lengths and reasons.
- ``dy47073_comparison.tsv``: per-gene three-way comparison (optional).
- ``dy47073_only.tsv``: DY47073 proteins outside the PomBase coding set (optional).
- ``afdb_only_dy47073.tsv``: AFdb-only entries matched to DY47073 (optional).

Usage
-----
    DY=tmp/afdb_vs_pombase/reference/dy47073
    python compare_afdb_pombase.py \
        --afdb-dir /hugedata/YushengYang/AFdb/UP000002485_284812_SCHPO_v6 \
        --pae-dir /hugedata/YushengYang/AFdb/PAE/UP000002485_284812_SCHPO_v6 \
        --pombase-dir tmp/afdb_vs_pombase/reference/pombase-2026-11-01 \
        --output-dir tmp/afdb_vs_pombase/outputs \
        --dy47073-fasta "$DY/S_pombe_DY47073.liftoff_braker3_polished_240717_protein.fasta" \
        --uniprot-resolution tmp/afdb_vs_pombase/intermediate/uniprot_resolution.tsv \
        --uniprot-triangulation tmp/afdb_vs_pombase/intermediate/uniprot_triangulation.tsv

Author: Yusheng Yang (guidance) + Agent (implementation)
Date:   2026-11-04
Version: 1.0.0
"""

# =============================================================================
# 1. IMPORTS (Standard Library -> Third-Party -> Local Project)
# =============================================================================
import argparse
from dataclasses import dataclass
from datetime import date
from enum import StrEnum
import gzip
from pathlib import Path
import re
import sys
from typing import TextIO

from Bio.Align import PairwiseAligner
from loguru import logger
import numpy as np
import pandas as pd

# =============================================================================
# 2. CONSTANTS & ENUMS (Domain Markers & Column Definitions)
# =============================================================================
PROTEIN_CODING = "protein coding gene"
ACCESSION_COLUMN_ALIASES = ("external_id", "uniprot_id")
CIF_SEQUENCE_KEY = "_entity_poly.pdbx_seq_one_letter_code_can"
CIF_DESCRIPTION_KEY = "_entity.pdbx_description"
CIF_NAME_PATTERN = re.compile(r"^AF-(?P<accession>.+)-F\d+-model_v\d+\.cif\.gz$")
PAE_NAME_PATTERN = re.compile(r"^AF-(?P<accession>.+)-F\d+-predicted_aligned_error_v\d+\.json$")
UNIPROT_RESOLUTION_COLUMNS = [
    "query_accession",
    "is_secondary",
    "primary_accession",
    "primary_in_pombase",
    "reviewed",
    "length",
    "protein_name",
]
UNIPROT_TRIANGULATION_COLUMNS = [
    "gene_systematic_id",
    "uniprot_id",
    "relationship",
    "uniprot_length",
]
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

SEQUENCE_CATEGORIES = (
    "identical",
    "substitution_only",
    "indel_only",
    "both",
)
AFDB_ONLY_STATUSES = (
    "sequence_identical_to_shared_protein",
    "represents_pombase_only_gene",
    "sequence_identical_to_other_afdb_entry",
    "no_pombase_match",
)
POMBASE_ONLY_COLUMNS = [
    "gene_systematic_id",
    "gene_name",
    "gene_product",
    "uniprot_id",
    "pombase_length",
    "matched_afdb_only_accessions",
    "reason",
    "exceeds_afdb_max_length",
]
MERGED_COLUMNS = [
    "category",
    "match_method",
    "uniprot_id",
    "afdb_accession",
    "pombase_accession",
    "gene_systematic_id",
    "gene_name",
    "gene_product",
    "gene_type",
    "afdb_description",
    "afdb_length",
    "pombase_length",
    "length_delta",
    "sequence_category",
    "n_substitutions",
    "n_afdb_extra",
    "n_pombase_extra",
    "n_afdb_extra_terminal",
    "n_pombase_extra_terminal",
    "identity_pct",
    "pae_available",
    "duplicate_afdb_accessions",
    "notes",
]
EXCEL_ENGINE = "openpyxl"
INTEGER_COLUMNS = [
    "afdb_length",
    "pombase_length",
    "length_delta",
    "n_substitutions",
    "n_afdb_extra",
    "n_pombase_extra",
    "n_afdb_extra_terminal",
    "n_pombase_extra_terminal",
]
DY47073_COMPARISON_COLUMNS = [
    "gene_systematic_id",
    "gene_name",
    "gene_product",
    "gene_type",
    "uniprot_id",
    "afdb_accession",
    "afdb_length",
    "pombase_length",
    "dy47073_length",
    "dy47073_n_records",
    "dy47073_record_ids",
    "dy47073_internal_stop_truncated",
    "dy47073_match_pombase",
    "dy47073_match_afdb",
    "dy47073_length_delta",
    "dy47073_vs_pombase",
    "dy47073_vs_afdb",
    "relationship",
    "dy47073_category",
    "dy47073_notes",
]
DY47073_ONLY_COLUMNS = [
    "gene_systematic_id",
    "gene_name",
    "gene_product",
    "gene_type",
    "dy47073_record_ids",
    "dy47073_length",
    "dy47073_internal_stop_truncated",
    "dy47073_matching_afdb_accessions",
    "notes",
]
AFDB_ONLY_DY_COLUMNS = [
    "uniprot_id",
    "afdb_description",
    "afdb_length",
    "dy47073_matching_genes",
    "dy47073_length",
    "dy47073_internal_stop_truncated",
    "dy47073_category",
]
DY47073_CATEGORY_LABELS = {
    "identical": "三者一致或 DY47073 == PomBase(无 AFdb)",
    "dy_eq_pombase_ne_afdb": "DY47073 == PomBase,AFdb 不同",
    "dy_eq_afdb_ne_pombase": "DY47073 == AFdb,PomBase 不同",
    "dy_differs_from_both": "DY47073 与两者都不同",
    "dy_differs_pombase_no_afdb": "DY47073 与 PomBase 不同,且无 AFdb 结构",
    "absent_from_dy47073": "DY47073 无该蛋白",
}


class GeneCol(StrEnum):
    """PomBase gene table column names."""

    SYSTEMATIC_ID = "gene_systematic_id"
    GENE_NAME = "gene_name"
    PRODUCT = "gene_product"
    GENE_TYPE = "gene_type"
    ACCESSION = "uniprot_id"


class AfdbCol(StrEnum):
    """AFdb structure table column names."""

    ACCESSION = "uniprot_id"
    DESCRIPTION = "afdb_description"
    SEQUENCE = "afdb_sequence"
    LENGTH = "afdb_length"


class PeptideCol(StrEnum):
    """PomBase peptide table column names."""

    SYSTEMATIC_ID = "gene_systematic_id"
    SEQUENCE = "pombase_sequence"
    LENGTH = "pombase_length"
    INTERNAL_STOP = "internal_stop_truncated"


class ComparisonCol(StrEnum):
    """Shared-protein comparison column names."""

    ACCESSION = "uniprot_id"
    AFDB_LENGTH = "afdb_length"
    POMBASE_LENGTH = "pombase_length"
    LENGTH_DELTA = "length_delta"
    CATEGORY = "sequence_category"
    N_SUBSTITUTIONS = "n_substitutions"
    N_AFDB_EXTRA = "n_afdb_extra"
    N_POMBASE_EXTRA = "n_pombase_extra"
    N_AFDB_EXTRA_TERMINAL = "n_afdb_extra_terminal"
    N_POMBASE_EXTRA_TERMINAL = "n_pombase_extra_terminal"
    IDENTITY_PCT = "identity_pct"


class DyCol(StrEnum):
    """DY47073 protein and three-way comparison column names."""

    GENE = "gene_systematic_id"
    RECORD_ID = "dy47073_record_id"
    SEQUENCE = "dy47073_sequence"
    LENGTH = "dy47073_length"
    N_RECORDS = "dy47073_n_records"
    RECORD_IDS = "dy47073_record_ids"
    INTERNAL_STOP = "dy47073_internal_stop_truncated"
    MATCH_POMBASE = "dy47073_match_pombase"
    MATCH_AFDB = "dy47073_match_afdb"
    LENGTH_DELTA = "dy47073_length_delta"
    VS_POMBASE = "dy47073_vs_pombase"
    VS_AFDB = "dy47073_vs_afdb"
    MATCHING_AFDB = "dy47073_matching_afdb_accessions"
    MATCHING_GENES = "dy47073_matching_genes"
    CATEGORY = "dy47073_category"
    NOTES = "dy47073_notes"


# =============================================================================
# 3. CONFIGURATION (Immutable Runtime Parameters)
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class Config:
    """Immutable runtime configuration container."""

    afdb_dir: Path
    pae_dir: Path
    pombase_dir: Path
    output_dir: Path
    pombase_release: str
    uniprot_resolution: Path | None = None
    uniprot_triangulation: Path | None = None
    dy47073_fasta: Path | None = None
    dy47073_label: str = "240717"


@dataclass(kw_only=True, frozen=True)
class Dy47073Result:
    """Container for the DY47073 three-way comparison tables."""

    comparison: pd.DataFrame
    only: pd.DataFrame
    afdb_only_matches: pd.DataFrame
    secondary_matches: pd.DataFrame
    record_count: int
    internal_stop_count: int


# =============================================================================
# 4. BOUNDARY LOADERS (I/O Shell)
# =============================================================================
def open_text(path: Path) -> TextIO:
    """Open a plain-text or gzip-compressed text file."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def read_cif_header_sequence(handle: TextIO) -> tuple[str, str]:
    """Read the entity description and modelled sequence from an AFdb CIF header."""
    description = ""
    sequence_chunks: list[str] = []
    in_sequence = False

    for line in handle:
        if in_sequence:
            if line.strip() == ";":
                break
            sequence_chunks.append(line.strip().lstrip(";"))
        elif line.startswith(CIF_SEQUENCE_KEY):
            remainder = line[len(CIF_SEQUENCE_KEY) :].strip()
            if remainder and not remainder.startswith(";"):
                sequence_chunks.append(remainder)
                break
            if remainder:
                sequence_chunks.append(remainder[1:])
            in_sequence = True
        elif line.startswith(CIF_DESCRIPTION_KEY):
            description = line.split(None, 1)[1].strip().strip('"')

    return description, "".join(sequence_chunks)


def parse_afdb_cif(path: Path) -> dict[str, str]:
    """Extract accession, protein description, and modelled sequence from one AFdb CIF."""
    name_match = CIF_NAME_PATTERN.match(path.name)
    if name_match is None:
        raise ValueError(f"Unrecognised AFdb CIF filename: {path.name}")

    with gzip.open(path, "rt") as handle:
        description, sequence = read_cif_header_sequence(handle)
    if not sequence:
        raise ValueError(f"Empty sequence parsed from {path.name}")

    return {
        AfdbCol.ACCESSION: name_match.group("accession"),
        AfdbCol.DESCRIPTION: description,
        AfdbCol.SEQUENCE: sequence,
    }


def load_afdb_entries(afdb_dir: Path) -> pd.DataFrame:
    """Load every AFdb model as one row of accession, description, and sequence."""
    paths = sorted(afdb_dir.glob("AF-*-F*-model_v*.cif.gz"))
    if not paths:
        raise FileNotFoundError(f"No AFdb CIF files found under {afdb_dir}")

    entries = pd.DataFrame.from_records([parse_afdb_cif(path) for path in paths])
    entries[AfdbCol.LENGTH] = entries[AfdbCol.SEQUENCE].str.len()
    unexpected = set("".join(entries[AfdbCol.SEQUENCE])) - STANDARD_AMINO_ACIDS
    if unexpected:
        logger.warning(f"Non-standard residues in AFdb sequences: {sorted(unexpected)}")
    logger.info(f"Parsed {len(entries)} AFdb CIF files from {afdb_dir}")
    return entries


def load_pae_accessions(pae_dir: Path) -> set[str]:
    """Return the accessions that have a PAE JSON file."""
    paths = sorted(pae_dir.glob("AF-*-F*-predicted_aligned_error_v*.json"))
    if not paths:
        raise FileNotFoundError(f"No PAE JSON files found under {pae_dir}")

    accessions: set[str] = set()
    for path in paths:
        name_match = PAE_NAME_PATTERN.match(path.name)
        if name_match is None:
            raise ValueError(f"Unrecognised PAE filename: {path.name}")
        accessions.add(name_match.group("accession"))
    logger.info(f"Found {len(accessions)} PAE JSON files under {pae_dir}")
    return accessions


def load_uniprot_resolution(path: Path, pombase_accessions: set[str]) -> pd.DataFrame:
    """Load the optional UniProt accession resolution table against a PomBase release."""
    resolution = pd.read_csv(path, sep="\t", dtype=str)
    missing = set(UNIPROT_RESOLUTION_COLUMNS) - set(resolution.columns)
    if missing:
        raise ValueError(f"Missing UniProt resolution columns: {sorted(missing)}")
    resolution = resolution.loc[:, UNIPROT_RESOLUTION_COLUMNS]
    resolution = resolution.rename(columns={"query_accession": AfdbCol.ACCESSION.value})
    resolution["primary_in_pombase"] = (
        resolution["primary_accession"].isin(pombase_accessions).astype(str)
    )
    return resolution


def load_uniprot_triangulation(path: Path) -> pd.DataFrame:
    """Load the optional current-UniProt triangulation table."""
    triangulation = pd.read_csv(path, sep="\t", dtype=str)
    missing = set(UNIPROT_TRIANGULATION_COLUMNS) - set(triangulation.columns)
    if missing:
        raise ValueError(f"Missing UniProt triangulation columns: {sorted(missing)}")
    triangulation = triangulation.loc[:, UNIPROT_TRIANGULATION_COLUMNS]
    triangulation["uniprot_length"] = pd.to_numeric(
        triangulation["uniprot_length"], errors="coerce"
    ).astype("Int64")
    return triangulation


def resolve_accession_column(columns: pd.Index) -> str:
    """Return the gene-table column that holds UniProt accessions."""
    for name in ACCESSION_COLUMN_ALIASES:
        if name in columns:
            return name
    raise KeyError(f"None of {ACCESSION_COLUMN_ALIASES} found in the gene table")


def load_pombase_gene_catalog(pombase_dir: Path) -> pd.DataFrame:
    """Load every PomBase gene, including non-coding genes, with its UniProt accession."""
    gene_table = pd.read_csv(pombase_dir / "gene_IDs_names_products.tsv", sep="\t", dtype=str)
    accession_column = resolve_accession_column(gene_table.columns)

    required = {GeneCol.SYSTEMATIC_ID, GeneCol.GENE_NAME, GeneCol.PRODUCT, GeneCol.GENE_TYPE}
    missing = required - set(gene_table.columns)
    if missing:
        raise ValueError(f"Missing gene-table columns: {sorted(missing)}")

    catalog = gene_table.rename(columns={accession_column: GeneCol.ACCESSION})
    return catalog.loc[:, list(GeneCol)].reset_index(drop=True)


def truncate_peptide_at_first_stop(raw_sequence: str) -> tuple[str, bool]:
    """Truncate a peptide at the first stop codon and report whether the stop was internal."""
    sequence, separator, tail = raw_sequence.partition("*")
    return sequence, bool(separator and tail)


def read_peptide_fasta(path: Path) -> list[tuple[str, str]]:
    """Read raw peptide records from a PomBase FASTA file."""
    raw_records: list[tuple[str, str]] = []
    systematic_id = ""
    chunks: list[str] = []

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                if systematic_id:
                    raw_records.append((systematic_id, "".join(chunks)))
                header = line[1:].split()[0].split(":pep")[0]
                systematic_id = header.rsplit(".", 1)[0]
                chunks = []
            else:
                chunks.append(line)
    if systematic_id:
        raw_records.append((systematic_id, "".join(chunks)))
    return raw_records


def load_pombase_peptides(pombase_dir: Path) -> pd.DataFrame:
    """Parse the PomBase peptide FASTA into gene-level sequences cut at the first stop codon."""
    fasta_path = pombase_dir / "peptide.fa.gz"
    if not fasta_path.exists():
        fasta_path = pombase_dir / "peptide.fa"

    raw_records = read_peptide_fasta(fasta_path)
    records = []
    for systematic_id, raw_sequence in raw_records:
        sequence, had_internal_stop = truncate_peptide_at_first_stop(raw_sequence)
        if had_internal_stop:
            logger.warning(
                f"Internal stop codon in PomBase peptide {systematic_id}; truncated at first stop"
            )
        records.append((systematic_id, sequence, had_internal_stop))

    peptides = pd.DataFrame(
        records,
        columns=pd.Index(
            [
                PeptideCol.SYSTEMATIC_ID.value,
                PeptideCol.SEQUENCE.value,
                PeptideCol.INTERNAL_STOP.value,
            ]
        ),
    )
    peptides[PeptideCol.LENGTH] = peptides[PeptideCol.SEQUENCE].str.len()

    is_duplicated = peptides[PeptideCol.SYSTEMATIC_ID].duplicated()
    if is_duplicated.any():
        duplicated_ids = peptides.loc[is_duplicated, PeptideCol.SYSTEMATIC_ID].tolist()
        raise ValueError(f"Multiple peptides for one gene: {duplicated_ids}")

    unexpected = set("".join(peptides[PeptideCol.SEQUENCE])) - STANDARD_AMINO_ACIDS
    if unexpected:
        logger.warning(f"Non-standard residues in PomBase peptides: {sorted(unexpected)}")
    return peptides


def read_dy47073_fasta(path: Path) -> list[tuple[str, str]]:
    """Read raw DY47073 protein records from a FASTA file."""
    raw_records: list[tuple[str, str]] = []
    record_id = ""
    chunks: list[str] = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                if record_id:
                    raw_records.append((record_id, "".join(chunks)))
                record_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if record_id:
        raw_records.append((record_id, "".join(chunks)))
    return raw_records


def resolve_dy47073_gene(record_id: str, gene_ids: set[str]) -> str | None:
    """Resolve one DY47073 record ID to a PomBase gene systematic ID."""
    header = record_id.removeprefix("braker_")
    candidates = (
        header,
        re.sub(r"_\d+\.\(1\+\)$", "", header),
        re.sub(r"_(?:1|2)$", "", header),
        re.sub(r"\.\d+$", "", re.sub(r"_(?:1|2)$", "", header)),
        re.sub(r"\.\d+$", "", re.sub(r"_\d+\.\(1\+\)$", "", header)),
    )
    return next((candidate for candidate in candidates if candidate in gene_ids), None)


def load_dy47073_proteins(path: Path, gene_catalog: pd.DataFrame) -> pd.DataFrame:
    """Parse DY47073 proteins into record-level rows with resolved PomBase gene IDs."""
    gene_ids = set(gene_catalog[GeneCol.SYSTEMATIC_ID])
    records = []
    for record_id, raw_sequence in read_dy47073_fasta(path):
        gene = resolve_dy47073_gene(record_id, gene_ids)
        records.append((gene or record_id, record_id, raw_sequence.rstrip("*")))

    proteins = pd.DataFrame(
        records,
        columns=pd.Index([DyCol.GENE.value, DyCol.RECORD_ID.value, DyCol.SEQUENCE.value]),
    )
    proteins[DyCol.LENGTH] = proteins[DyCol.SEQUENCE].str.len()
    proteins[DyCol.INTERNAL_STOP] = proteins[DyCol.SEQUENCE].str.contains("*", regex=False)

    unexpected = set("".join(proteins[DyCol.SEQUENCE])) - STANDARD_AMINO_ACIDS - {"*"}
    if unexpected:
        logger.warning(f"Non-standard residues in DY47073 proteins: {sorted(unexpected)}")
    internal_stop_count = int(proteins[DyCol.INTERNAL_STOP].sum())
    if internal_stop_count:
        logger.warning(f"{internal_stop_count} DY47073 records contain internal stop codons")
    logger.info(f"Parsed {len(proteins)} DY47073 protein records from {path}")
    return proteins


# =============================================================================
# 5. FUNCTIONAL CORE (Pure, Trusted Transformations)
# =============================================================================
def build_aligner() -> PairwiseAligner:
    """Build a global protein aligner for difference classification."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def gap_position_counts(target: str, query: str) -> tuple[int, int]:
    """Count target residues opposite query gaps at terminal versus internal positions."""
    gap_positions = [
        index
        for index, (target_residue, query_residue) in enumerate(zip(target, query, strict=True))
        if target_residue != "-" and query_residue == "-"
    ]
    paired_positions = [
        index
        for index, (target_residue, query_residue) in enumerate(zip(target, query, strict=True))
        if target_residue != "-" and query_residue != "-"
    ]
    if not gap_positions or not paired_positions:
        return 0, len(gap_positions)

    first_paired = paired_positions[0]
    last_paired = paired_positions[-1]
    terminal = sum(1 for index in gap_positions if index < first_paired or index > last_paired)
    return terminal, len(gap_positions) - terminal


def compare_sequences(
    afdb_sequence: str,
    pombase_sequence: str,
    aligner: PairwiseAligner,
) -> dict[str, int | float | str]:
    """Classify one protein pair as identical, substitution-only, indel-only, or both."""
    if afdb_sequence == pombase_sequence:
        return {
            ComparisonCol.CATEGORY: "identical",
            ComparisonCol.N_SUBSTITUTIONS: 0,
            ComparisonCol.N_AFDB_EXTRA: 0,
            ComparisonCol.N_POMBASE_EXTRA: 0,
            ComparisonCol.N_AFDB_EXTRA_TERMINAL: 0,
            ComparisonCol.N_POMBASE_EXTRA_TERMINAL: 0,
            ComparisonCol.IDENTITY_PCT: 100.0,
        }

    alignment = aligner.align(afdb_sequence, pombase_sequence)[0]
    afdb_aligned = alignment[0]
    pombase_aligned = alignment[1]

    n_substitutions = sum(
        1
        for afdb_residue, pombase_residue in zip(afdb_aligned, pombase_aligned, strict=True)
        if afdb_residue != pombase_residue and "-" not in (afdb_residue, pombase_residue)
    )
    n_afdb_extra = pombase_aligned.count("-")
    n_pombase_extra = afdb_aligned.count("-")
    n_matches = sum(
        1
        for afdb_residue, pombase_residue in zip(afdb_aligned, pombase_aligned, strict=True)
        if afdb_residue == pombase_residue and afdb_residue != "-"
    )

    if n_afdb_extra == 0 and n_pombase_extra == 0:
        category = "substitution_only"
    elif n_substitutions == 0:
        category = "indel_only"
    else:
        category = "both"

    afdb_extra_terminal, _ = gap_position_counts(afdb_aligned, pombase_aligned)
    pombase_extra_terminal, _ = gap_position_counts(pombase_aligned, afdb_aligned)
    paired_residues = n_matches + n_substitutions
    identity_pct = round(100.0 * n_matches / paired_residues, 2) if paired_residues else 0.0

    return {
        ComparisonCol.CATEGORY: category,
        ComparisonCol.N_SUBSTITUTIONS: n_substitutions,
        ComparisonCol.N_AFDB_EXTRA: n_afdb_extra,
        ComparisonCol.N_POMBASE_EXTRA: n_pombase_extra,
        ComparisonCol.N_AFDB_EXTRA_TERMINAL: afdb_extra_terminal,
        ComparisonCol.N_POMBASE_EXTRA_TERMINAL: pombase_extra_terminal,
        ComparisonCol.IDENTITY_PCT: identity_pct,
    }


def split_accession_sets(
    afdb_entries: pd.DataFrame,
    pombase_genes: pd.DataFrame,
) -> tuple[set[str], set[str], set[str]]:
    """Split accessions into shared, AFdb-only, and PomBase-only sets."""
    afdb_accessions = set(afdb_entries[AfdbCol.ACCESSION])
    pombase_accessions = set(pombase_genes[GeneCol.ACCESSION].dropna())
    return (
        afdb_accessions & pombase_accessions,
        afdb_accessions - pombase_accessions,
        pombase_accessions - afdb_accessions,
    )


def attach_comparison(
    shared_genes: pd.DataFrame,
    afdb_entries: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    aligner: PairwiseAligner,
) -> pd.DataFrame:
    """Attach lengths and sequence-difference statistics to each shared gene/accession pair."""
    merged = shared_genes.merge(afdb_entries, on=AfdbCol.ACCESSION, how="inner")
    merged = merged.merge(pombase_peptides, on=PeptideCol.SYSTEMATIC_ID, how="inner")
    merged[ComparisonCol.LENGTH_DELTA] = merged[AfdbCol.LENGTH] - merged[PeptideCol.LENGTH]

    statistics = pd.DataFrame.from_records(
        [
            compare_sequences(afdb_sequence, pombase_sequence, aligner)
            for afdb_sequence, pombase_sequence in zip(
                merged[AfdbCol.SEQUENCE], merged[PeptideCol.SEQUENCE], strict=True
            )
        ]
    )
    return pd.concat([merged.reset_index(drop=True), statistics], axis=1)


def join_values(mapping: dict[str, list[str]]) -> dict[str, str]:
    """Join sequence-keyed list values into semicolon-separated strings."""
    return {sequence: ";".join(values) for sequence, values in mapping.items()}


def annotate_afdb_only(
    afdb_only: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    shared_entries: pd.DataFrame,
    shared_genes: pd.DataFrame,
) -> pd.DataFrame:
    """Match AFdb-only entries to PomBase peptides and shared AFdb sequences by exact sequence."""
    peptide_to_genes = (
        pombase_peptides.groupby(PeptideCol.SEQUENCE)[PeptideCol.SYSTEMATIC_ID].agg(list).to_dict()
    )
    shared_sequence_to_accessions = (
        shared_entries.groupby(AfdbCol.SEQUENCE)[AfdbCol.ACCESSION].agg(list).to_dict()
    )
    shared_gene_ids = set(shared_genes[GeneCol.SYSTEMATIC_ID])
    peptide_to_shared = {
        sequence: [gene for gene in genes if gene in shared_gene_ids]
        for sequence, genes in peptide_to_genes.items()
    }
    peptide_to_pombase_only = {
        sequence: [gene for gene in genes if gene not in shared_gene_ids]
        for sequence, genes in peptide_to_genes.items()
    }

    annotated = afdb_only.assign(
        matched_pombase_genes=afdb_only[AfdbCol.SEQUENCE]
        .map(join_values(peptide_to_genes))
        .fillna(""),
        matched_shared_genes=afdb_only[AfdbCol.SEQUENCE]
        .map(join_values(peptide_to_shared))
        .fillna(""),
        matched_pombase_only_genes=afdb_only[AfdbCol.SEQUENCE]
        .map(join_values(peptide_to_pombase_only))
        .fillna(""),
        duplicate_of_shared_afdb_accessions=afdb_only[AfdbCol.SEQUENCE]
        .map(join_values(shared_sequence_to_accessions))
        .fillna(""),
    )
    status_conditions = [
        annotated["matched_shared_genes"] != "",
        annotated["matched_pombase_only_genes"] != "",
        annotated["duplicate_of_shared_afdb_accessions"] != "",
    ]
    annotated["match_status"] = np.select(
        status_conditions,
        [
            "sequence_identical_to_shared_protein",
            "represents_pombase_only_gene",
            "sequence_identical_to_other_afdb_entry",
        ],
        default="no_pombase_match",
    )
    return annotated


def build_pombase_only(
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    afdb_only: pd.DataFrame,
    afdb_max_length: int,
) -> pd.DataFrame:
    """List PomBase genes that have no direct AFdb accession match."""
    genes_with_peptides = pombase_genes.merge(
        pombase_peptides, on=PeptideCol.SYSTEMATIC_ID, how="inner"
    )
    joined_matches = join_values(
        afdb_only.groupby(AfdbCol.SEQUENCE)[AfdbCol.ACCESSION].agg(list).to_dict()
    )
    genes_without_accession = genes_with_peptides.loc[
        genes_with_peptides[GeneCol.ACCESSION].isna()
    ].copy()
    genes_without_accession["uniprot_id"] = ""
    genes_without_accession["matched_afdb_only_accessions"] = (
        genes_without_accession[PeptideCol.SEQUENCE].map(joined_matches).fillna("")
    )
    genes_without_accession["reason"] = "no_uniprot_accession"
    genes_without_accession["exceeds_afdb_max_length"] = (
        genes_without_accession[PeptideCol.LENGTH] > afdb_max_length
    )

    is_represented = genes_without_accession["matched_afdb_only_accessions"] != ""
    genes_without_accession.loc[is_represented, "reason"] = (
        "represented_in_afdb_under_other_accession"
    )
    return genes_without_accession.loc[:, POMBASE_ONLY_COLUMNS]


def build_pombase_only_missing_accessions(
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    pombase_only_accessions: set[str],
    afdb_only: pd.DataFrame,
    afdb_max_length: int,
) -> pd.DataFrame:
    """List PomBase-only accessions with length and exact-sequence reconciliation."""
    joined_matches = join_values(
        afdb_only.groupby(AfdbCol.SEQUENCE)[AfdbCol.ACCESSION].agg(list).to_dict()
    )
    genes = pombase_genes.loc[pombase_genes[GeneCol.ACCESSION].isin(pombase_only_accessions)]
    genes = genes.merge(pombase_peptides, on=PeptideCol.SYSTEMATIC_ID, how="inner")
    genes["matched_afdb_only_accessions"] = (
        genes[PeptideCol.SEQUENCE].map(joined_matches).fillna("")
    )
    genes["reason"] = "not_modelled_by_afdb"
    genes["exceeds_afdb_max_length"] = genes[PeptideCol.LENGTH] > afdb_max_length
    genes.loc[genes["exceeds_afdb_max_length"], "reason"] = "exceeds_afdb_length_cutoff"

    is_represented = genes["matched_afdb_only_accessions"] != ""
    genes.loc[is_represented, "reason"] = "represented_in_afdb_under_other_accession"
    return genes.loc[:, POMBASE_ONLY_COLUMNS]


def uniprot_triangulation_counts(shared_comparison: pd.DataFrame) -> dict[str, int]:
    """Summarize the optional AFdb/PomBase/current-UniProt triangulation."""
    if "relationship" not in shared_comparison.columns:
        return {}

    relationship = shared_comparison["relationship"]
    return {
        "uniprot_check_pombase_matches": int(
            (relationship == "pombase_matches_uniprot_afdb_differs").sum()
        ),
        "uniprot_check_afdb_matches": int(
            (relationship == "afdb_matches_uniprot_pombase_differs").sum()
        ),
        "uniprot_check_neither": int((relationship == "neither_matches_uniprot").sum()),
        "uniprot_check_missing": int((relationship == "uniprot_missing").sum()),
    }


def uniprot_resolution_counts(afdb_only: pd.DataFrame) -> dict[str, int]:
    """Summarise the optional UniProt accession resolution."""
    if "is_secondary" not in afdb_only.columns:
        return {}

    is_secondary = afdb_only["is_secondary"] == "True"
    primary_in_pombase = afdb_only["primary_in_pombase"] == "True"
    return {
        "uniprot_resolved": int(afdb_only["is_secondary"].notna().sum()),
        "uniprot_secondary": int(is_secondary.sum()),
        "uniprot_secondary_in_pombase": int((is_secondary & primary_in_pombase).sum()),
        "uniprot_primary": int((~is_secondary).sum()),
        "uniprot_primary_in_pombase": int(((~is_secondary) & primary_in_pombase).sum()),
    }


def build_summary_counts(
    afdb_entries: pd.DataFrame,
    pae_accessions: set[str],
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    shared_accessions: set[str],
    afdb_only: pd.DataFrame,
    pombase_only_accessions: set[str],
    shared_comparison: pd.DataFrame,
    pombase_only: pd.DataFrame,
) -> dict[str, int]:
    """Collect headline counts for the report."""
    categories = shared_comparison[ComparisonCol.CATEGORY]
    afdb_accessions = set(afdb_entries[AfdbCol.ACCESSION])
    counts = {
        "afdb_entries": len(afdb_entries),
        "afdb_unique_accessions": afdb_entries[AfdbCol.ACCESSION].nunique(),
        "pombase_coding_genes": len(pombase_genes),
        "pombase_genes_with_accession": int(pombase_genes[GeneCol.ACCESSION].notna().sum()),
        "pombase_unique_accessions": int(pombase_genes[GeneCol.ACCESSION].nunique()),
        "shared_accessions": len(shared_accessions),
        "afdb_only_accessions": len(afdb_only),
        "pombase_only_accessions": len(pombase_only_accessions),
        "pombase_genes_without_accession": int(pombase_genes[GeneCol.ACCESSION].isna().sum()),
        "shared_genes": len(shared_comparison),
        "pombase_only_genes": len(pombase_only),
        "same_length_pairs": int(
            (
                shared_comparison[ComparisonCol.AFDB_LENGTH]
                == shared_comparison[ComparisonCol.POMBASE_LENGTH]
            ).sum()
        ),
        "length_diff_pairs": int(
            (
                shared_comparison[ComparisonCol.AFDB_LENGTH]
                != shared_comparison[ComparisonCol.POMBASE_LENGTH]
            ).sum()
        ),
        "identical_pairs": int((categories == "identical").sum()),
        "substitution_only_pairs": int((categories == "substitution_only").sum()),
        "indel_only_pairs": int((categories == "indel_only").sum()),
        "both_pairs": int((categories == "both").sum()),
        "afdb_only_matching_shared": int(
            (afdb_only["match_status"] == "sequence_identical_to_shared_protein").sum()
        ),
        "afdb_only_representing_pombase_only": int(
            (afdb_only["match_status"] == "represents_pombase_only_gene").sum()
        ),
        "afdb_only_no_match": int((afdb_only["match_status"] == "no_pombase_match").sum()),
        "afdb_max_length": int(afdb_entries[AfdbCol.LENGTH].max()),
        "pombase_max_length": int(pombase_peptides[PeptideCol.LENGTH].max()),
        "pombase_max_shared_length": int(shared_comparison[ComparisonCol.POMBASE_LENGTH].max()),
        "pae_entries": len(pae_accessions),
        "afdb_accessions_without_pae": len(afdb_accessions - pae_accessions),
        "pae_accessions_without_afdb": len(pae_accessions - afdb_accessions),
        "pombase_internal_stop_peptides": int(pombase_peptides[PeptideCol.INTERNAL_STOP].sum()),
        "pombase_only_represented": int(
            (pombase_only["reason"] == "represented_in_afdb_under_other_accession").sum()
        ),
        "pombase_only_length_cutoff": int(
            (pombase_only["reason"] == "exceeds_afdb_length_cutoff").sum()
        ),
        "pombase_only_absent_other": int(
            pombase_only["reason"].isin(["not_modelled_by_afdb", "no_uniprot_accession"]).sum()
        ),
    }
    counts.update(uniprot_resolution_counts(afdb_only))
    return counts


def join_unique(values: pd.Series) -> str:
    """Join the unique non-missing string values of a Series with semicolons."""
    return ";".join(sorted({str(value) for value in values if pd.notna(value)}))


def build_shared_by_sequence(
    afdb_only: pd.DataFrame,
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
) -> pd.DataFrame:
    """Build shared rows for PomBase genes matched to AFdb entries by exact sequence."""
    matched = afdb_only.loc[afdb_only["match_status"] == "represents_pombase_only_gene"].copy()
    pairs = matched.assign(
        gene_systematic_id=matched["matched_pombase_only_genes"].str.split(";")
    ).explode("gene_systematic_id")
    pairs = pairs.loc[pairs["gene_systematic_id"] != ""]

    grouped = (
        pairs.groupby("gene_systematic_id")
        .agg(
            afdb_accession=(AfdbCol.ACCESSION.value, join_unique),
            afdb_description=(AfdbCol.DESCRIPTION.value, "first"),
            afdb_length=(AfdbCol.LENGTH.value, "first"),
        )
        .reset_index()
    )

    genes = pombase_genes.rename(columns={GeneCol.ACCESSION: "pombase_accession"})
    peptides = pombase_peptides.loc[:, [PeptideCol.SYSTEMATIC_ID, PeptideCol.LENGTH]]
    grouped = grouped.merge(genes, on=GeneCol.SYSTEMATIC_ID, how="inner")
    grouped = grouped.merge(peptides, on=PeptideCol.SYSTEMATIC_ID, how="inner")

    grouped["category"] = "shared"
    grouped["match_method"] = "sequence"
    grouped["uniprot_id"] = grouped["afdb_accession"].str.split(";").str[0]
    grouped["length_delta"] = 0
    grouped["sequence_category"] = "identical"
    grouped["n_substitutions"] = 0
    grouped["n_afdb_extra"] = 0
    grouped["n_pombase_extra"] = 0
    grouped["n_afdb_extra_terminal"] = 0
    grouped["n_pombase_extra_terminal"] = 0
    grouped["identity_pct"] = 100.0
    grouped["pae_available"] = True
    grouped["duplicate_afdb_accessions"] = ""
    grouped["notes"] = "matched_by_sequence; pombase_accession_not_in_afdb"
    has_no_accession = grouped["pombase_accession"].isna()
    grouped.loc[has_no_accession, "notes"] = "matched_by_sequence; pombase_has_no_uniprot_accession"
    return grouped


def add_shared_match_columns(
    shared_comparison: pd.DataFrame,
    afdb_only: pd.DataFrame,
) -> pd.DataFrame:
    """Add match provenance and duplicate AFdb accessions to ID-matched shared rows."""
    duplicates = afdb_only.loc[
        afdb_only["match_status"] == "sequence_identical_to_shared_protein"
    ].copy()
    duplicate_pairs = duplicates.assign(
        gene_systematic_id=duplicates["matched_shared_genes"].str.split(";")
    ).explode("gene_systematic_id")
    duplicate_pairs = duplicate_pairs.loc[duplicate_pairs["gene_systematic_id"] != ""]
    duplicate_map = (
        duplicate_pairs.groupby("gene_systematic_id")[AfdbCol.ACCESSION.value]
        .agg(join_unique)
        .to_dict()
    )

    annotated = shared_comparison.assign(
        category="shared",
        match_method="uniprot_id",
        afdb_accession=shared_comparison[ComparisonCol.ACCESSION],
        pombase_accession=shared_comparison[ComparisonCol.ACCESSION],
        pae_available=True,
        duplicate_afdb_accessions=shared_comparison[GeneCol.SYSTEMATIC_ID]
        .map(duplicate_map)
        .fillna(""),
    )
    annotated["notes"] = ""
    has_duplicates = annotated["duplicate_afdb_accessions"] != ""
    annotated.loc[has_duplicates, "notes"] = "afdb_secondary_accession"
    return annotated


def build_afdb_only_rows(afdb_only: pd.DataFrame) -> pd.DataFrame:
    """Build main-table rows for AFdb entries with no PomBase protein match."""
    unmatched = afdb_only.loc[afdb_only["match_status"] == "no_pombase_match"].copy()
    notes = pd.Series("no_pombase_match", index=unmatched.index)
    if "is_secondary" in unmatched.columns:
        is_secondary = unmatched["is_secondary"] == "True"
        notes = notes.mask(is_secondary, "no_pombase_match; uniprot_secondary")

    return unmatched.assign(
        category="afdb_only",
        match_method="none",
        uniprot_id=unmatched[AfdbCol.ACCESSION],
        afdb_accession=unmatched[AfdbCol.ACCESSION],
        pombase_accession=np.nan,
        gene_systematic_id=np.nan,
        gene_name=np.nan,
        gene_product=np.nan,
        gene_type=np.nan,
        pombase_length=np.nan,
        length_delta=np.nan,
        sequence_category="not_compared",
        n_substitutions=np.nan,
        n_afdb_extra=np.nan,
        n_pombase_extra=np.nan,
        n_afdb_extra_terminal=np.nan,
        n_pombase_extra_terminal=np.nan,
        identity_pct=np.nan,
        pae_available=True,
        duplicate_afdb_accessions="",
        notes=notes,
    )


def build_pombase_only_rows(pombase_only: pd.DataFrame) -> pd.DataFrame:
    """Build main-table rows for PomBase genes with no AFdb structure."""
    unrepresented = pombase_only.loc[
        pombase_only["reason"] != "represented_in_afdb_under_other_accession"
    ].copy()
    return unrepresented.assign(
        category="pombase_only",
        match_method="none",
        uniprot_id=unrepresented["uniprot_id"],
        afdb_accession=np.nan,
        pombase_accession=unrepresented["uniprot_id"],
        gene_type=PROTEIN_CODING,
        afdb_description=np.nan,
        afdb_length=np.nan,
        length_delta=np.nan,
        sequence_category="not_compared",
        n_substitutions=np.nan,
        n_afdb_extra=np.nan,
        n_pombase_extra=np.nan,
        n_afdb_extra_terminal=np.nan,
        n_pombase_extra_terminal=np.nan,
        identity_pct=np.nan,
        pae_available=False,
        duplicate_afdb_accessions="",
        notes=unrepresented["reason"],
    )


def build_merged_table(
    shared_by_id: pd.DataFrame,
    shared_by_sequence: pd.DataFrame,
    afdb_only_rows: pd.DataFrame,
    pombase_only_rows: pd.DataFrame,
) -> pd.DataFrame:
    """Concatenate the unified protein table and select the report columns."""
    merged = pd.concat(
        [shared_by_id, shared_by_sequence, afdb_only_rows, pombase_only_rows],
        ignore_index=True,
    ).reindex(columns=MERGED_COLUMNS)
    for column in INTEGER_COLUMNS:
        merged[column] = merged[column].astype("Int64")
    return merged


def build_afdb_sequence_table(
    shared_by_id: pd.DataFrame,
    afdb_only: pd.DataFrame,
) -> pd.DataFrame:
    """Build one AFdb sequence per gene, preferring the primary accession."""
    matched = afdb_only.loc[
        afdb_only["match_status"].isin(
            ["sequence_identical_to_shared_protein", "represents_pombase_only_gene"]
        )
    ]
    secondary_rows = []
    for gene_column in ("matched_shared_genes", "matched_pombase_only_genes"):
        pairs = matched.assign(gene_systematic_id=matched[gene_column].str.split(";")).explode(
            "gene_systematic_id"
        )
        is_gene = pairs["gene_systematic_id"].notna() & (pairs["gene_systematic_id"] != "")
        secondary_rows.append(
            pairs.loc[
                is_gene,
                ["gene_systematic_id", AfdbCol.ACCESSION, AfdbCol.SEQUENCE, AfdbCol.LENGTH],
            ].rename(columns={AfdbCol.ACCESSION: "afdb_accession"})
        )

    primary_rows = shared_by_id.loc[
        :, ["gene_systematic_id", "afdb_accession", AfdbCol.SEQUENCE, AfdbCol.LENGTH]
    ]
    combined = pd.concat([*secondary_rows, primary_rows], ignore_index=True)
    # Primary rows come last, so they win for genes that also have a secondary entry.
    return combined.drop_duplicates("gene_systematic_id", keep="last").reset_index(drop=True)


def build_dy47073_gene_table(
    dy47073_proteins: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    afdb_sequences: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate DY47073 records per gene and flag exact matches to PomBase and AFdb."""
    annotated = dy47073_proteins.merge(
        pombase_peptides.loc[:, [PeptideCol.SYSTEMATIC_ID, PeptideCol.SEQUENCE, PeptideCol.LENGTH]],
        on=GeneCol.SYSTEMATIC_ID,
        how="left",
    )
    annotated = annotated.merge(afdb_sequences, on=GeneCol.SYSTEMATIC_ID, how="left")
    annotated["match_pombase"] = annotated[DyCol.SEQUENCE] == annotated[PeptideCol.SEQUENCE]
    annotated["match_afdb"] = annotated[DyCol.SEQUENCE] == annotated[AfdbCol.SEQUENCE]

    grouped = annotated.groupby(DyCol.GENE, as_index=False).agg(
        **{
            DyCol.N_RECORDS: (DyCol.RECORD_ID, "size"),
            DyCol.RECORD_IDS: (DyCol.RECORD_ID, join_unique),
            DyCol.LENGTH: (DyCol.LENGTH, "max"),
            DyCol.INTERNAL_STOP: (DyCol.INTERNAL_STOP, "any"),
            DyCol.MATCH_POMBASE: ("match_pombase", "any"),
            DyCol.MATCH_AFDB: ("match_afdb", "any"),
        }
    )
    representative = (
        dy47073_proteins.sort_values(DyCol.LENGTH, ascending=False, kind="stable")
        .drop_duplicates(DyCol.GENE)
        .loc[:, [DyCol.GENE, DyCol.SEQUENCE]]
    )
    return grouped.merge(representative, on=DyCol.GENE, how="left")


def build_dy47073_notes(comparison: pd.DataFrame) -> pd.Series:
    """Summarise DY47073 isoform and internal-stop flags for one row per gene."""
    notes = pd.Series("", index=comparison.index, dtype="object")
    multiple = comparison[DyCol.N_RECORDS] > 1
    internal_stop = comparison[DyCol.INTERNAL_STOP]
    notes = notes.mask(internal_stop, "internal_stop")
    has_notes = notes.ne("")
    notes = notes.mask(multiple & ~has_notes, "multiple_records")
    return notes.mask(multiple & has_notes, notes + ";multiple_records")


def build_dy47073_comparison(
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    afdb_sequences: pd.DataFrame,
    dy47073_genes: pd.DataFrame,
) -> pd.DataFrame:
    """Build the per-gene three-way comparison for every protein-coding PomBase gene."""
    base = pombase_genes.merge(
        pombase_peptides.loc[:, [PeptideCol.SYSTEMATIC_ID, PeptideCol.LENGTH]],
        on=GeneCol.SYSTEMATIC_ID,
    )
    base = base.merge(afdb_sequences, on=GeneCol.SYSTEMATIC_ID, how="left")
    base = base.rename(columns={PeptideCol.LENGTH: "pombase_length"})

    dy_columns = [
        GeneCol.SYSTEMATIC_ID,
        DyCol.N_RECORDS,
        DyCol.RECORD_IDS,
        DyCol.LENGTH,
        DyCol.SEQUENCE,
        DyCol.INTERNAL_STOP,
        DyCol.MATCH_POMBASE,
        DyCol.MATCH_AFDB,
    ]
    comparison = base.merge(dy47073_genes.loc[:, dy_columns], on=GeneCol.SYSTEMATIC_ID, how="left")
    comparison[DyCol.N_RECORDS] = comparison[DyCol.N_RECORDS].fillna(0).astype("Int64")
    for column in (DyCol.INTERNAL_STOP, DyCol.MATCH_POMBASE, DyCol.MATCH_AFDB):
        comparison[column] = comparison[column].astype("boolean").fillna(False).astype(bool)
    comparison[DyCol.RECORD_IDS] = comparison[DyCol.RECORD_IDS].fillna("")

    has_dy = comparison[DyCol.N_RECORDS] > 0
    has_afdb = comparison[AfdbCol.SEQUENCE].notna()
    matches_pombase = comparison[DyCol.MATCH_POMBASE]
    matches_afdb = comparison[DyCol.MATCH_AFDB]

    comparison[DyCol.VS_POMBASE] = np.select(
        [~has_dy, matches_pombase, comparison[DyCol.LENGTH] == comparison["pombase_length"]],
        ["absent_from_dy47073", "identical", "same_length_different"],
        default="different_length",
    )
    comparison[DyCol.VS_AFDB] = np.select(
        [
            ~has_afdb,
            ~has_dy,
            matches_afdb,
            comparison[DyCol.LENGTH] == comparison[AfdbCol.LENGTH],
        ],
        ["no_afdb_structure", "absent_from_dy47073", "identical", "same_length_different"],
        default="different_length",
    )
    comparison[DyCol.CATEGORY] = np.select(
        [
            ~has_dy,
            matches_pombase & (~has_afdb | matches_afdb),
            matches_pombase & has_afdb & ~matches_afdb,
            ~matches_pombase & matches_afdb,
            ~matches_pombase & ~has_afdb,
        ],
        [
            "absent_from_dy47073",
            "identical",
            "dy_eq_pombase_ne_afdb",
            "dy_eq_afdb_ne_pombase",
            "dy_differs_pombase_no_afdb",
        ],
        default="dy_differs_from_both",
    )
    comparison[DyCol.LENGTH_DELTA] = comparison[DyCol.LENGTH] - comparison["pombase_length"]
    for column in (AfdbCol.LENGTH, DyCol.LENGTH, DyCol.LENGTH_DELTA):
        comparison[column] = comparison[column].astype("Int64")
    comparison[DyCol.NOTES] = build_dy47073_notes(comparison)
    return comparison


def build_dy47073_only_table(
    dy47073_genes: pd.DataFrame,
    dy47073_proteins: pd.DataFrame,
    gene_catalog: pd.DataFrame,
    afdb_only: pd.DataFrame,
) -> pd.DataFrame:
    """Build the table of DY47073 proteins outside the PomBase protein-coding set."""
    coding_ids = set(
        gene_catalog.loc[gene_catalog[GeneCol.GENE_TYPE] == PROTEIN_CODING, GeneCol.SYSTEMATIC_ID]
    )
    only = dy47073_genes.loc[~dy47073_genes[DyCol.GENE].isin(coding_ids)].copy()
    metadata = gene_catalog.set_index(GeneCol.SYSTEMATIC_ID)[
        [GeneCol.GENE_NAME, GeneCol.PRODUCT, GeneCol.GENE_TYPE]
    ]
    only = only.merge(metadata, left_on=DyCol.GENE, right_index=True, how="left")
    only[GeneCol.GENE_TYPE] = only[GeneCol.GENE_TYPE].fillna("not_in_pombase_2026")

    dy_sequences = dy47073_proteins.loc[:, [DyCol.GENE, DyCol.SEQUENCE]].rename(
        columns={DyCol.SEQUENCE: AfdbCol.SEQUENCE}
    )
    matches = dy_sequences.merge(
        afdb_only.loc[:, [AfdbCol.SEQUENCE, AfdbCol.ACCESSION]],
        on=AfdbCol.SEQUENCE,
        how="inner",
    )
    accessions = matches.groupby(DyCol.GENE)[AfdbCol.ACCESSION].agg(join_unique)
    only[DyCol.MATCHING_AFDB] = only[DyCol.GENE].map(accessions).fillna("")

    only["notes"] = "pombase_2026_pseudogene"
    only.loc[only[GeneCol.GENE_TYPE] == "not_in_pombase_2026", "notes"] = "not_in_pombase_2026"
    return only


def annotate_afdb_only_with_dy47073(
    afdb_only: pd.DataFrame,
    dy47073_proteins: pd.DataFrame,
) -> pd.DataFrame:
    """Match AFdb-only entries without a PomBase hit to DY47073 proteins by exact sequence."""
    unmatched = afdb_only.loc[afdb_only["match_status"] == "no_pombase_match"].copy()
    dy_sequences = dy47073_proteins.loc[
        :, [DyCol.GENE, DyCol.SEQUENCE, DyCol.LENGTH, DyCol.INTERNAL_STOP]
    ].rename(columns={DyCol.SEQUENCE: AfdbCol.SEQUENCE})
    matches = dy_sequences.merge(
        unmatched.loc[:, [AfdbCol.SEQUENCE, AfdbCol.ACCESSION]],
        on=AfdbCol.SEQUENCE,
        how="inner",
    )
    grouped = matches.groupby(AfdbCol.ACCESSION).agg(
        **{
            DyCol.MATCHING_GENES: (DyCol.GENE, join_unique),
            DyCol.LENGTH: (DyCol.LENGTH, "max"),
            DyCol.INTERNAL_STOP: (DyCol.INTERNAL_STOP, "any"),
        }
    )
    annotated = unmatched.merge(grouped, on=AfdbCol.ACCESSION, how="left")
    annotated[DyCol.MATCHING_GENES] = annotated[DyCol.MATCHING_GENES].fillna("")
    annotated[DyCol.INTERNAL_STOP] = (
        annotated[DyCol.INTERNAL_STOP].astype("boolean").fillna(False).astype(bool)
    )
    annotated[DyCol.LENGTH] = annotated[DyCol.LENGTH].astype("Int64")
    annotated[DyCol.CATEGORY] = np.where(
        annotated[DyCol.MATCHING_GENES] != "", "dy47073_sequence_match", "no_dy47073_match"
    )
    return annotated


def build_dy47073_secondary_matches(
    dy47073_proteins: pd.DataFrame,
    afdb_only: pd.DataFrame,
    dy_comparison: pd.DataFrame,
) -> pd.DataFrame:
    """Find DY47073 proteins that match an AFdb-only accession other than the gene's primary."""
    dy_sequences = dy47073_proteins.loc[
        :, [DyCol.GENE, DyCol.RECORD_ID, DyCol.SEQUENCE, DyCol.LENGTH]
    ].rename(columns={DyCol.SEQUENCE: AfdbCol.SEQUENCE})
    matches = dy_sequences.merge(
        afdb_only.loc[:, [AfdbCol.SEQUENCE, AfdbCol.ACCESSION, AfdbCol.LENGTH]],
        on=AfdbCol.SEQUENCE,
        how="inner",
    )
    primary = dy_comparison.set_index(GeneCol.SYSTEMATIC_ID)[["afdb_accession", DyCol.VS_POMBASE]]
    matches = matches.merge(primary, left_on=DyCol.GENE, right_index=True, how="left")
    is_coding = matches[DyCol.GENE].isin(dy_comparison[GeneCol.SYSTEMATIC_ID])
    is_secondary = matches["afdb_accession"].isna() | (
        matches[AfdbCol.ACCESSION] != matches["afdb_accession"]
    )
    differs_from_pombase = matches[DyCol.VS_POMBASE] != "identical"
    keep = is_coding & is_secondary & differs_from_pombase
    return matches.loc[keep].reset_index(drop=True)


def attach_dy47073_to_merged(
    merged: pd.DataFrame,
    dy_result: Dy47073Result,
) -> pd.DataFrame:
    """Add DY47073 status columns to the unified merged table."""
    dy_columns = [
        GeneCol.SYSTEMATIC_ID,
        DyCol.LENGTH,
        DyCol.N_RECORDS,
        DyCol.VS_POMBASE,
        DyCol.VS_AFDB,
        DyCol.CATEGORY,
    ]
    enriched = merged.merge(
        dy_result.comparison.loc[:, dy_columns], on=GeneCol.SYSTEMATIC_ID, how="left"
    )
    afdb_only_matches = dy_result.afdb_only_matches.set_index(AfdbCol.ACCESSION)
    is_afdb_only = enriched["category"] == "afdb_only"
    for column in (DyCol.CATEGORY, DyCol.LENGTH):
        enriched.loc[is_afdb_only, column] = enriched.loc[is_afdb_only, AfdbCol.ACCESSION].map(
            afdb_only_matches[column]
        )
    enriched[DyCol.CATEGORY] = enriched[DyCol.CATEGORY].fillna("not_compared")
    return enriched


def dy47073_counts(dy_result: Dy47073Result) -> dict[str, int]:
    """Collect headline counts for the DY47073 three-way comparison."""
    comparison = dy_result.comparison
    category = comparison[DyCol.CATEGORY]
    differs = category == "dy_differs_from_both"
    same_length = comparison[DyCol.LENGTH] == comparison[PeptideCol.LENGTH]
    return {
        "dy47073_records": dy_result.record_count,
        "dy47073_internal_stop_records": dy_result.internal_stop_count,
        "dy47073_coding_genes_with_record": int(comparison[DyCol.N_RECORDS].gt(0).sum()),
        "dy47073_identical": int((category == "identical").sum()),
        "dy47073_pombase_matches_afdb_differs": int((category == "dy_eq_pombase_ne_afdb").sum()),
        "dy47073_afdb_matches_pombase_differs": int((category == "dy_eq_afdb_ne_pombase").sum()),
        "dy47073_differs_from_both": int(differs.sum()),
        "dy47073_differs_same_length": int((differs & same_length).sum()),
        "dy47073_differs_different_length": int((differs & ~same_length).sum()),
        "dy47073_differs_from_pombase_no_afdb": int(
            (category == "dy_differs_pombase_no_afdb").sum()
        ),
        "dy47073_absent_from_pombase": int((category == "absent_from_dy47073").sum()),
        "dy47073_only_genes": len(dy_result.only),
        "dy47073_only_with_afdb": int((dy_result.only[DyCol.MATCHING_AFDB] != "").sum()),
        "afdb_only_dy47073_match": int(
            (dy_result.afdb_only_matches[DyCol.CATEGORY] == "dy47073_sequence_match").sum()
        ),
        "afdb_only_no_dy47073_match": int(
            (dy_result.afdb_only_matches[DyCol.CATEGORY] == "no_dy47073_match").sum()
        ),
    }


def build_afdb_entry_table(
    afdb_entries: pd.DataFrame,
    pae_accessions: set[str],
    shared_by_id: pd.DataFrame,
    shared_by_sequence: pd.DataFrame,
    afdb_only: pd.DataFrame,
) -> pd.DataFrame:
    """Map every AFdb entry to its protein entity and entry role."""
    roles = pd.concat(
        [
            build_primary_entry_roles(shared_by_id),
            build_duplicate_entry_roles(shared_by_id),
            build_sequence_entry_roles(shared_by_sequence),
            build_unmatched_entry_roles(afdb_only),
        ],
        ignore_index=True,
    )
    roles = roles.groupby("afdb_accession", as_index=False).agg(
        gene_systematic_id=("gene_systematic_id", join_unique),
        gene_name=("gene_name", join_unique),
        pombase_accession=("pombase_accession", join_unique),
        entry_role=("entry_role", "first"),
        protein_entity_accession=("protein_entity_accession", join_unique),
        notes=("notes", join_unique),
    )
    roles["entity_category"] = np.where(
        roles["entry_role"] == "no_pombase_match", "afdb_only", "shared"
    )

    roles = roles.rename(columns={"afdb_accession": AfdbCol.ACCESSION.value})
    table = afdb_entries.merge(roles, on=AfdbCol.ACCESSION, how="left")
    if table["entry_role"].isna().any():
        missing = table.loc[table["entry_role"].isna(), AfdbCol.ACCESSION].tolist()
        raise ValueError(f"AFdb entries without a protein-entity mapping: {missing[:10]}")

    table["pae_available"] = table[AfdbCol.ACCESSION].isin(pae_accessions)
    table = table.rename(columns={AfdbCol.ACCESSION: "afdb_accession"})
    return table.loc[
        :,
        [
            "afdb_accession",
            AfdbCol.DESCRIPTION.value,
            AfdbCol.LENGTH.value,
            "pae_available",
            "entity_category",
            "entry_role",
            "protein_entity_accession",
            "gene_systematic_id",
            "gene_name",
            "pombase_accession",
            "notes",
        ],
    ]


def build_primary_entry_roles(shared_by_id: pd.DataFrame) -> pd.DataFrame:
    """Build entry roles for the current AFdb primary accession of shared genes."""
    return shared_by_id.assign(
        entry_role="primary",
        protein_entity_accession=shared_by_id["afdb_accession"],
        notes="primary",
    ).loc[
        :,
        [
            "afdb_accession",
            "gene_systematic_id",
            "gene_name",
            "pombase_accession",
            "entry_role",
            "protein_entity_accession",
            "notes",
        ],
    ]


def build_duplicate_entry_roles(shared_by_id: pd.DataFrame) -> pd.DataFrame:
    """Build entry roles for AFdb duplicate entries listed on shared rows."""
    duplicates = shared_by_id.loc[shared_by_id["duplicate_afdb_accessions"] != ""].copy()
    duplicates["entry_role"] = "duplicate_secondary_accession"
    duplicates["protein_entity_accession"] = duplicates["afdb_accession"]
    duplicates["notes"] = "duplicate of " + duplicates["protein_entity_accession"]
    duplicates = duplicates.assign(
        afdb_accession=duplicates["duplicate_afdb_accessions"].str.split(";")
    ).explode("afdb_accession")
    return duplicates.loc[
        :,
        [
            "afdb_accession",
            "gene_systematic_id",
            "gene_name",
            "pombase_accession",
            "entry_role",
            "protein_entity_accession",
            "notes",
        ],
    ]


def build_sequence_entry_roles(shared_by_sequence: pd.DataFrame) -> pd.DataFrame:
    """Build entry roles for AFdb entries matched to PomBase genes by sequence."""
    sequence = shared_by_sequence.copy()
    sequence["entry_role"] = "sequence_matched"
    sequence["protein_entity_accession"] = sequence["uniprot_id"]
    sequence["notes"] = "matched_by_sequence"
    sequence = sequence.assign(afdb_accession=sequence["afdb_accession"].str.split(";")).explode(
        "afdb_accession"
    )
    return sequence.loc[
        :,
        [
            "afdb_accession",
            "gene_systematic_id",
            "gene_name",
            "pombase_accession",
            "entry_role",
            "protein_entity_accession",
            "notes",
        ],
    ]


def build_unmatched_entry_roles(afdb_only: pd.DataFrame) -> pd.DataFrame:
    """Build entry roles for AFdb entries with no PomBase match."""
    unmatched = afdb_only.loc[afdb_only["match_status"] == "no_pombase_match"].copy()
    unmatched["entry_role"] = "no_pombase_match"
    unmatched["protein_entity_accession"] = unmatched[AfdbCol.ACCESSION]
    unmatched["gene_systematic_id"] = ""
    unmatched["gene_name"] = ""
    unmatched["pombase_accession"] = ""
    unmatched["notes"] = "no_pombase_match"
    if "is_secondary" in unmatched.columns:
        is_secondary = unmatched["is_secondary"] == "True"
        unmatched.loc[is_secondary, "notes"] = "no_pombase_match; uniprot_secondary"
    return unmatched.loc[
        :,
        [
            AfdbCol.ACCESSION,
            "gene_systematic_id",
            "gene_name",
            "pombase_accession",
            "entry_role",
            "protein_entity_accession",
            "notes",
        ],
    ].rename(columns={AfdbCol.ACCESSION: "afdb_accession"})


# =============================================================================
# 6. REPORT RENDERING
# =============================================================================
def format_table_cell(value: object) -> str:
    """Format one Markdown table cell, replacing missing values and escaping pipes."""
    if pd.isna(value):
        return ""
    return str(value).replace("|", "\\|").replace("\n", " ")


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    """Render a compact Markdown table from selected DataFrame columns."""
    header = "| " + " | ".join(columns) + " |"
    separator = "|" + "|".join("---" for _ in columns) + "|"
    body = [
        "| " + " | ".join(format_table_cell(value) for value in row) + " |"
        for row in frame.loc[:, columns].itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *body])


def render_view(
    frame: pd.DataFrame,
    columns: list[str],
    headers: list[str],
    sort_column: str | None = None,
    limit: int | None = None,
) -> str:
    """Render an optionally sorted and truncated DataFrame view as a Markdown table."""
    view = frame.sort_values(sort_column, ascending=False) if sort_column else frame
    if limit is not None:
        view = view.head(limit)
    view = view.loc[:, columns]
    view.columns = headers
    return markdown_table(view, headers)


def render_set_overview(counts: dict[str, int], config: Config) -> str:
    """Render the set-membership overview section."""
    rows = [
        ("AFdb 结构条目(全部 F1 片段)", counts["afdb_entries"]),
        ("AFdb 唯一 UniProt accession", counts["afdb_unique_accessions"]),
        ("PAE JSON 文件", counts["pae_entries"]),
        ("PomBase protein coding genes", counts["pombase_coding_genes"]),
        ("PomBase 有 UniProt accession 的基因", counts["pombase_genes_with_accession"]),
        ("PomBase 唯一 UniProt accession", counts["pombase_unique_accessions"]),
        ("共有 accession", counts["shared_accessions"]),
        ("AFdb 特有 accession(含重复/序列匹配)", counts["afdb_only_accessions"]),
        ("PomBase 特有 accession", counts["pombase_only_accessions"]),
        ("PomBase 无 UniProt accession 的基因", counts["pombase_genes_without_accession"]),
        ("共有基因(UniProt ID 匹配)", counts["shared_genes"]),
        ("共有基因(序列匹配)", counts["shared_by_sequence_genes"]),
        ("共有基因合计(合并表 shared)", counts["merged_shared_rows"]),
        ("AFdb 特有蛋白(合并表 afdb_only)", counts["afdb_only_proteins"]),
        ("PomBase 特有基因(合并表 pombase_only)", counts["pombase_only_proteins"]),
        ("合并表总行数", counts["merged_rows"]),
    ]
    table = "\n".join(f"| {name} | {value:,} |" for name, value in rows)
    return f"""## 1. 集合总览

数据版本:

- AFdb:`UP000002485_284812_SCHPO_v6`,{counts["afdb_entries"]:,} 个 CIF 结构文件。
- PomBase:`{config.pombase_release}`,{counts["pombase_coding_genes"]:,} 个 protein coding gene。
- AFdb 序列取 CIF 的 `_entity_poly.pdbx_seq_one_letter_code_can`(实际建模序列)。
- PomBase 序列取 `peptide.fa`(在第一个终止符处截断;本版本有
  {counts["pombase_internal_stop_peptides"]} 条内部终止符记录)。
- 主匹配键为 UniProt accession;PomBase 无 accession 或 accession 不在 AFdb 时,
  再用 peptide 与 AFdb 序列做精确匹配。合并表 `match_method` 列记录 `uniprot_id` 或 `sequence`。

| 项目 | 数量 |
|---|---:|
{table}
"""


def render_length_section(shared_comparison: pd.DataFrame, counts: dict[str, int]) -> str:
    """Render the length-agreement section."""
    length_diff = shared_comparison.loc[
        shared_comparison[ComparisonCol.AFDB_LENGTH]
        != shared_comparison[ComparisonCol.POMBASE_LENGTH]
    ].copy()
    length_diff["abs_delta"] = length_diff[ComparisonCol.LENGTH_DELTA].abs()
    length_diff = length_diff.sort_values("abs_delta", ascending=False)
    top = length_diff.head(15).loc[
        :,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            ComparisonCol.ACCESSION,
            ComparisonCol.AFDB_LENGTH,
            ComparisonCol.POMBASE_LENGTH,
            ComparisonCol.LENGTH_DELTA,
        ],
    ]
    top.columns = ["gene", "name", "UniProt", "AFdb 长度", "PomBase 长度", "差值"]
    table = markdown_table(top, list(top.columns))
    return f"""## 2. 共有蛋白的长度一致性

以下统计基于 {counts["shared_genes"]:,} 个 UniProt ID 匹配的共有基因对;
另有 {counts["shared_by_sequence_genes"]:,} 个序列匹配的共有基因对,均为等长且完全一致。

- 长度相同:{counts["same_length_pairs"]:,} / {counts["shared_genes"]:,} 个共有基因对。
- 长度不同:{counts["length_diff_pairs"]:,} / {counts["shared_genes"]:,} 个共有基因对。
- AFdb 最长建模序列:{counts["afdb_max_length"]:,} aa。
- PomBase 全蛋白组最长:{counts["pombase_max_length"]:,} aa;
  共有蛋白中 PomBase 最长:{counts["pombase_max_shared_length"]:,} aa。

长度差异最大的 15 个共有基因对:

{table}

注:最大的几例(gpi14、sap34、rhp16)经 UniProt 复核,AFdb 序列与当前 UniProt 长度一致;
差异来自 UniProt 与 PomBase 基因模型不一致,而不是 AFdb 截断。
"""


def render_sequence_section(shared_comparison: pd.DataFrame, counts: dict[str, int]) -> str:
    """Render the sequence-agreement section."""
    identical = counts["identical_pairs"]
    substitution = counts["substitution_only_pairs"]
    indel = counts["indel_only_pairs"]
    both = counts["both_pairs"]

    nonidentical = shared_comparison.loc[
        shared_comparison[ComparisonCol.CATEGORY] != "identical"
    ].copy()
    nonidentical["difference_score"] = (
        nonidentical[ComparisonCol.N_SUBSTITUTIONS]
        + nonidentical[ComparisonCol.N_AFDB_EXTRA]
        + nonidentical[ComparisonCol.N_POMBASE_EXTRA]
    )
    nonidentical = nonidentical.sort_values("difference_score", ascending=False)
    top = nonidentical.head(30).loc[
        :,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            ComparisonCol.ACCESSION,
            ComparisonCol.AFDB_LENGTH,
            ComparisonCol.POMBASE_LENGTH,
            ComparisonCol.CATEGORY,
            ComparisonCol.N_SUBSTITUTIONS,
            ComparisonCol.N_AFDB_EXTRA,
            ComparisonCol.N_POMBASE_EXTRA,
            ComparisonCol.IDENTITY_PCT,
        ],
    ]
    top.columns = [
        "gene",
        "name",
        "UniProt",
        "AFdb",
        "PomBase",
        "类型",
        "替换",
        "AFdb 多出",
        "AFdb 缺少",
        "identity %",
    ]
    table = markdown_table(top, list(top.columns))

    substitutions = shared_comparison.loc[
        shared_comparison[ComparisonCol.CATEGORY] == "substitution_only"
    ].sort_values(ComparisonCol.N_SUBSTITUTIONS, ascending=False)
    substitution_view = substitutions.loc[
        :,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            ComparisonCol.ACCESSION,
            ComparisonCol.N_SUBSTITUTIONS,
            ComparisonCol.IDENTITY_PCT,
        ],
    ]
    substitution_view.columns = ["gene", "name", "UniProt", "替换数", "identity %"]
    substitution_table = markdown_table(substitution_view, list(substitution_view.columns))

    if "uniprot_check_pombase_matches" in counts:
        triangulation_note = f"""三方核对(2026-11-04 当前 UniProt):

- {counts["uniprot_check_pombase_matches"]:,} 个:PomBase == 当前 UniProt,AFdb 不同
  → AFdb v6 使用了旧序列,当前参考应以 PomBase/UniProt 为准。
- {counts["uniprot_check_afdb_matches"]:,} 个:AFdb == 当前 UniProt,PomBase 不同
  → UniProt 与 PomBase 注释分歧,需逐条核对基因模型。
- {counts["uniprot_check_neither"]:,} 个:三者都不一致,需人工核对。
"""
    else:
        triangulation_note = (
            "本次额外对照了 2023-07-01、2024-03-01、2025-05-01、2025-06-01、2026-11-01 "
            "五个 PomBase release:102 个非一致蛋白中有 86 个的 AFdb 长度不等于任何一个 "
            "release 的 peptide 长度,说明差异主要不是 PomBase 版本更新造成的。\n"
        )

    return f"""## 3. 共有蛋白的序列一致性

以下统计基于 {counts["shared_genes"]:,} 个 UniProt ID 匹配的共有基因对;
另有 {counts["shared_by_sequence_genes"]:,} 个序列匹配对均为 identical,
加上后完全一致共 {identical + counts["shared_by_sequence_genes"]:,} 个。

| 类型 | 数量 |
|---|---:|
| identical(完全一致) | {identical:,} |
| substitution_only(等长,仅点突变) | {substitution:,} |
| indel_only(仅序列增减) | {indel:,} |
| both(点突变 + 序列增减) | {both:,} |

等长、仅点突变的蛋白({substitution:,} 个):

{substitution_table}

非一致蛋白共 {len(nonidentical):,} 个。下表按差异位点数排序,
列出前 30 个(完整表见 `nonidentical_pairs.tsv`):

{table}

说明:`identity %` 按可比对残基计算(不计 gap)。PomBase peptide 中的内部终止符已在第一个终止符处截断
({counts["pombase_internal_stop_peptides"]} 条),这些基因的"缺失"残基不应视为真实 indel。

{triangulation_note}"""


def render_afdb_only_section(afdb_only: pd.DataFrame, counts: dict[str, int]) -> str:
    """Render the AFdb-only reconciliation section."""
    status_counts = afdb_only["match_status"].value_counts()
    status_rows = "\n".join(
        f"| {status} | {int(status_counts.get(status, 0)):,} |" for status in AFDB_ONLY_STATUSES
    )
    unmatched = afdb_only.loc[afdb_only["match_status"] == "no_pombase_match"].copy()
    unmatched = unmatched.sort_values(AfdbCol.LENGTH, ascending=False)
    unmatched_view = unmatched.loc[
        :, [AfdbCol.ACCESSION, AfdbCol.DESCRIPTION, AfdbCol.LENGTH]
    ].head(30)
    unmatched_view.columns = ["UniProt", "描述", "长度"]
    unmatched_table = markdown_table(unmatched_view, list(unmatched_view.columns))

    matched = afdb_only.loc[afdb_only["match_status"] != "no_pombase_match"].copy()
    matched_view = matched.loc[
        :,
        [
            AfdbCol.ACCESSION,
            AfdbCol.LENGTH,
            "match_status",
            "matched_shared_genes",
            "matched_pombase_only_genes",
        ],
    ].head(30)
    matched_view.columns = ["UniProt", "长度", "状态", "匹配的共有基因", "匹配的 PomBase-only 基因"]
    matched_table = markdown_table(matched_view, list(matched_view.columns))

    uniprot_note = ""
    if "uniprot_resolved" in counts:
        stale_secondary = (
            counts["uniprot_secondary_in_pombase"] - counts["afdb_only_matching_shared"]
        )
        uniprot_note = (
            "UniProt accession 复核(独立于序列匹配):\n\n"
            f"- {counts['uniprot_resolved']:,} 个 AFdb-only accession "
            "全部解析到当前 UniProt 条目。\n"
            f"- {counts['uniprot_secondary']:,} 个是次级 accession;"
            f"{counts['uniprot_secondary_in_pombase']:,} 个的当前 primary 在 PomBase 中,"
            f"其中 {counts['afdb_only_matching_shared']:,} 个与当前 peptide 完全一致,"
            f"{stale_secondary:,} 个序列版本较旧。\n"
            f"- {counts['uniprot_primary']:,} 个是当前 primary 条目,"
            "其 accession 均不在 PomBase 的 protein-coding 映射中。\n"
        )

    return f"""## 4. AFdb 特有蛋白(accession 层面)

AFdb-only accession 共 {counts["afdb_only_accessions"]:,} 个。
按 AFdb 序列与 PomBase peptide 的精确匹配分类:

| 状态 | 数量 |
|---|---:|
{status_rows}

{uniprot_note}
- `sequence_identical_to_shared_protein`:该 AFdb 条目的序列与某个共有基因的 peptide 完全相同,
  通常说明 AFdb 使用了旧/次级 UniProt accession,或该蛋白存在等序列重复条目。
- `represents_pombase_only_gene`:该 AFdb 条目对应一个 accession 未被 AFdb 使用的 PomBase 基因,
  说明 PomBase-only 列表中的这个基因其实有结构。
- `no_pombase_match`:序列与所有 PomBase peptide 都不完全相同,
  是真正需要单独核对的 AFdb-only 条目。

序列匹配上的前 30 个 AFdb-only 条目:

{matched_table}

无 PomBase 精确匹配的 AFdb-only 条目(按长度降序,前 30 个):

{unmatched_table}
"""


def render_pombase_only_section(pombase_only: pd.DataFrame, counts: dict[str, int]) -> str:
    """Render the PomBase-only section."""
    view = pombase_only.sort_values("pombase_length", ascending=False).loc[
        :,
        [
            "gene_systematic_id",
            "gene_name",
            "gene_product",
            "uniprot_id",
            "pombase_length",
            "reason",
            "exceeds_afdb_max_length",
        ],
    ]
    view.columns = ["gene", "name", "product", "UniProt", "长度", "原因", "超过 AFdb 最大长度"]
    table = markdown_table(view, list(view.columns))
    return f"""## 5. PomBase 特有蛋白

PomBase 侧共有 {counts["pombase_only_genes"]:,} 个基因没有直接的 AFdb accession 匹配
(其中 {counts["pombase_only_accessions"]:,} 个 accession 不在 AFdb,
{counts["pombase_genes_without_accession"]:,} 个基因无 UniProt accession)。

- 其中 {counts["pombase_only_represented"]:,} 个在 AFdb 中以另一个 accession 存在
  (见 `matched_afdb_only_accessions` 列),已作为 `match_method=sequence` 计入合并表的
  shared 类别;下表仍按 accession 层面列出,便于核对。
- {counts["pombase_only_length_cutoff"]:,} 个超过 AFdb 最大建模长度
  ({counts["afdb_max_length"]:,} aa),AFdb 不建模。
- 其余 {counts["pombase_only_absent_other"]:,} 个低于长度上限但仍无 AFdb 结构,
  例如 B5BP45-49(在 UP000002485 中但 AFdb 未收录)和 mat2-Pc/Pi(不在 UP000002485 中)。

{table}
"""


def render_dy47073_section(
    dy_result: Dy47073Result,
    counts: dict[str, int],
    config: Config,
) -> str:
    """Render the DY47073 three-way comparison section."""
    comparison = dy_result.comparison
    category_counts = (
        comparison[DyCol.CATEGORY]
        .value_counts()
        .reindex(list(DY47073_CATEGORY_LABELS), fill_value=0)
    )
    category_table = "\n".join(
        f"| {DY47073_CATEGORY_LABELS[name]} | {int(category_counts[name]):,} |"
        for name in DY47073_CATEGORY_LABELS
    )
    overview_table = "\n".join(
        f"| {label} | {value:,} |"
        for label, value in (
            ("DY47073 蛋白记录", counts["dy47073_records"]),
            ("有内部终止符的记录", counts["dy47073_internal_stop_records"]),
            ("PomBase coding 基因中有 DY47073 记录", counts["dy47073_coding_genes_with_record"]),
            ("DY47073 特有基因(PomBase 非 coding)", counts["dy47073_only_genes"]),
            ("其中在 AFdb 中有结构", counts["dy47073_only_with_afdb"]),
            ("AFdb-only 条目中与 DY47073 精确匹配", counts["afdb_only_dy47073_match"]),
        )
    )

    pombase_matches = comparison.loc[comparison[DyCol.CATEGORY] == "dy_eq_pombase_ne_afdb"].assign(
        abs_afdb_delta=lambda frame: (frame[AfdbCol.LENGTH] - frame["pombase_length"]).abs()
    )
    pombase_table = render_view(
        pombase_matches,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            GeneCol.ACCESSION,
            "pombase_length",
            AfdbCol.LENGTH,
            DyCol.LENGTH,
            "abs_afdb_delta",
        ],
        ["gene", "name", "UniProt", "PomBase", "AFdb", "DY47073", "abs(AFdb-PomBase)"],
        sort_column="abs_afdb_delta",
        limit=15,
    )
    pombase_uniprot = int(
        (pombase_matches["relationship"] == "pombase_matches_uniprot_afdb_differs").sum()
    )
    afdb_uniprot = int(
        (pombase_matches["relationship"] == "afdb_matches_uniprot_pombase_differs").sum()
    )

    afdb_matches = comparison.loc[comparison[DyCol.CATEGORY] == "dy_eq_afdb_ne_pombase"]
    afdb_table = render_view(
        afdb_matches,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            GeneCol.ACCESSION,
            "pombase_length",
            AfdbCol.LENGTH,
            DyCol.LENGTH,
            "relationship",
        ],
        ["gene", "name", "UniProt", "PomBase", "AFdb", "DY47073", "UniProt 关系"],
    )

    differs = comparison.loc[comparison[DyCol.CATEGORY] == "dy_differs_from_both"].assign(
        abs_dy_delta=lambda frame: (frame[DyCol.LENGTH] - frame["pombase_length"]).abs()
    )
    differs_table = render_view(
        differs,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            GeneCol.ACCESSION,
            "pombase_length",
            AfdbCol.LENGTH,
            DyCol.LENGTH,
            DyCol.VS_POMBASE,
            "abs_dy_delta",
        ],
        [
            "gene",
            "name",
            "UniProt",
            "PomBase",
            "AFdb",
            "DY47073",
            "DY vs PomBase",
            "abs(DY-PomBase)",
        ],
        sort_column="abs_dy_delta",
        limit=15,
    )
    differs_afdb_identical = int(differs["relationship"].isna().sum())
    differs_afdb_other = int(differs["relationship"].notna().sum())

    no_afdb = comparison.loc[comparison[DyCol.CATEGORY] == "dy_differs_pombase_no_afdb"]
    no_afdb_table = render_view(
        no_afdb,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            GeneCol.ACCESSION,
            "pombase_length",
            DyCol.LENGTH,
            DyCol.VS_POMBASE,
        ],
        ["gene", "name", "UniProt", "PomBase", "DY47073", "DY vs PomBase"],
    )

    absent = comparison.loc[comparison[DyCol.CATEGORY] == "absent_from_dy47073"]
    absent_table = render_view(
        absent,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            GeneCol.ACCESSION,
            "pombase_length",
            AfdbCol.LENGTH,
        ],
        ["gene", "name", "UniProt", "PomBase", "AFdb"],
    )

    only_table = render_view(
        dy_result.only,
        [
            GeneCol.SYSTEMATIC_ID,
            GeneCol.GENE_NAME,
            GeneCol.PRODUCT,
            GeneCol.GENE_TYPE,
            DyCol.RECORD_IDS,
            DyCol.LENGTH,
            DyCol.INTERNAL_STOP,
            DyCol.MATCHING_AFDB,
        ],
        [
            "gene",
            "name",
            "product",
            "PomBase 类型",
            "DY records",
            "DY 长度",
            "内部终止符",
            "匹配的 AFdb 条目",
        ],
        sort_column=DyCol.LENGTH,
    )

    matched = dy_result.afdb_only_matches.loc[
        dy_result.afdb_only_matches[DyCol.CATEGORY] == "dy47073_sequence_match"
    ]
    matched_table = render_view(
        matched,
        [
            AfdbCol.ACCESSION,
            AfdbCol.DESCRIPTION,
            AfdbCol.LENGTH,
            DyCol.MATCHING_GENES,
            DyCol.LENGTH,
        ],
        ["UniProt", "描述", "AFdb 长度", "DY47073 基因", "DY 长度"],
        sort_column=AfdbCol.LENGTH,
    )

    if dy_result.secondary_matches.empty:
        secondary_block = "未发现 DY47073 与 AFdb-only 次级条目的第三种序列版本。"
    else:
        secondary_table = render_view(
            dy_result.secondary_matches,
            [
                DyCol.GENE,
                DyCol.RECORD_ID,
                DyCol.LENGTH,
                AfdbCol.ACCESSION,
                AfdbCol.LENGTH,
                DyCol.VS_POMBASE,
            ],
            [
                "gene",
                "DY record",
                "DY 长度",
                "匹配的 AFdb-only 条目",
                "AFdb 长度",
                "DY vs PomBase",
            ],
            sort_column=DyCol.LENGTH,
        )
        secondary_block = f"""以下 DY47073 蛋白与一个非 primary 的 AFdb-only 条目精确匹配,
且与 PomBase 不同,构成第三种序列版本:

{secondary_table}"""

    return f"""## 6. DY47073 三方比较

数据:`{config.dy47073_fasta}`(release `{config.dy47073_label}`)。
DY47073 序列按 FASTA 翻译序列读取,仅去掉末尾终止符;内部终止符保留并在表中标记。
同一基因有多条记录时取最长序列作为代表,`dy47073_n_records` 记录条数。

| 项目 | 数量 |
|---|---:|
{overview_table}

| 三方关系 | 数量 |
|---|---:|
{category_table}

- `identical`:DY47073 与 PomBase 精确一致;若该基因有 primary AFdb 结构,也与 AFdb 一致。
- `dy_eq_pombase_ne_afdb`:DY47073 支持 PomBase,AFdb 序列不同(通常为旧序列)。
- `dy_eq_afdb_ne_pombase`:DY47073 支持 AFdb,PomBase 序列不同。
- `dy_differs_from_both`:DY47073 与 PomBase、AFdb 都不同。
- `dy_differs_pombase_no_afdb`:DY47073 与 PomBase 不同,且没有 AFdb 结构。
- `absent_from_dy47073`:PomBase 有该基因,DY47073 无对应蛋白。

### 6.1 DY47073 == PomBase,AFdb 不同({counts["dy47073_pombase_matches_afdb_differs"]})

这些基因支持 PomBase,AFdb 是旧序列或基因模型不同。
其中 {pombase_uniprot:,} 个当前 UniProt 也等于 PomBase,
{afdb_uniprot:,} 个当前 UniProt 等于 AFdb。
按 abs(AFdb-PomBase) 长度差排序的前 15 个:

{pombase_table}

### 6.2 DY47073 == AFdb,PomBase 不同({counts["dy47073_afdb_matches_pombase_differs"]})

{afdb_table}

注意:这 {counts["dy47073_afdb_matches_pombase_differs"]:,} 个中只有部分表示"PomBase 更新";
若当前 UniProt 等于 PomBase,则 DY47073 与 AFdb 保留的是旧序列。

### 6.3 DY47073 与两者都不同({counts["dy47073_differs_from_both"]})

- {differs_afdb_identical:,} 个 AFdb == PomBase,DY47073 是唯一差异(株特异或注释差异)。
- {differs_afdb_other:,} 个 AFdb 与 PomBase 本来就不同,构成三方分歧。
- 其中等长 {counts["dy47073_differs_same_length"]:,} 个、
  长度不同 {counts["dy47073_differs_different_length"]:,} 个(相对 PomBase)。

按 abs(DY-PomBase) 长度差排序的前 15 个:

{differs_table}

{secondary_block}

### 6.4 DY47073 与 PomBase 不同且无 AFdb 结构({counts["dy47073_differs_from_pombase_no_afdb"]})

{no_afdb_table}

### 6.5 DY47073 缺少的 PomBase 基因({counts["dy47073_absent_from_pombase"]})

{absent_table}

### 6.6 DY47073 特有蛋白 / PomBase pseudogene({counts["dy47073_only_genes"]})

这些 ID 在 PomBase 2026 中被标为 pseudogene(或不在基因表中),但 DY47073 有蛋白记录。
其中 {counts["dy47073_only_with_afdb"]:,} 个能在 AFdb 中找到精确匹配的结构。

{only_table}

### 6.7 AFdb-only 条目与 DY47073

{counts["afdb_only_no_match"]:,} 个无 PomBase 匹配的 AFdb-only 条目中,
{counts["afdb_only_dy47073_match"]:,} 个与 DY47073 精确匹配;
其余 {counts["afdb_only_no_dy47073_match"]:,} 个仍无匹配。

{matched_table}
"""


def render_recommendations(
    counts: dict[str, int],
    dy_result: Dy47073Result | None,
) -> str:
    """Render downstream-use recommendations."""
    nonidentical = (
        counts["substitution_only_pairs"] + counts["indel_only_pairs"] + counts["both_pairs"]
    )
    section_number = 7 if dy_result is not None else 6
    dy_block = ""
    if dy_result is not None:
        dy_block = f"""6. **DY47073 三方比较**:
   {counts["dy47073_pombase_matches_afdb_differs"]:,} 个基因 DY47073 == PomBase 但 AFdb 不同,
   {counts["dy47073_afdb_matches_pombase_differs"]:,} 个 DY47073 == AFdb 但 PomBase 不同,
   {counts["dy47073_differs_from_both"]:,} 个与两者都不同。
   若下游以 DY47073 为参考,需按 `dy47073_category` 区分坐标体系;现有 vDulab2 自定义预测
   基于旧 240322 注释,建议按本表更新。
"""
    return f"""## {section_number}. 对下游使用的建议

1. **主分析集**:使用 {counts["shared_accessions"]:,} 个 ID 匹配的共有 accession,
   对应 {counts["shared_genes"]:,} 个基因;加上 {counts["shared_by_sequence_genes"]:,} 个
   序列匹配基因,合并表共有 {counts["merged_shared_rows"]:,} 个 shared 行。
   下游以 `results/registry/gene_structure_registry.tsv` 的 `in_analysis` 为准。
2. **AFdb-only 条目**:{counts["afdb_only_matching_shared"]:,} 个与共有蛋白序列完全相同。
   其中 12 个与 primary 序列相同(真正的重复);62 个是序列与 PomBase 精确一致的替代条目,
   在注册表中作为这些基因的首选结构。{counts["afdb_only_no_match"]:,} 个无 PomBase
   精确匹配的条目不进入分析。
3. **序列不一致的共有蛋白**:共 {nonidentical:,} 个。注册表按"结构序列 == PomBase peptide"
   选择:62 个使用 AFdb 替代条目,38 个使用 vDulab,SPCC569.03 因预测序列与 PomBase 不符而排除。
   由于选择条件是序列精确一致,下游不再需要 UniProt/PomBase 坐标转换。
4. **PomBase-only 基因**:{counts["pombase_only_length_cutoff"]:,} 个长蛋白超过 AFdb 最大建模长度
   ({counts["afdb_max_length"]:,} aa)。其中 6 个有 vDulab、3 个有 AF3;另外 mat2-Pc/Pi 也有
   vDulab。sib1 只有 1-702 片段、B5BP45-49 无结构,均排除出分析。
5. **PAE 配对**:AFdb 结构 {counts["afdb_entries"]:,} 个,PAE 文件 {counts["pae_entries"]:,} 个;
   结构无 PAE {counts["afdb_accessions_without_pae"]:,} 个,PAE 无结构
   {counts["pae_accessions_without_afdb"]:,} 个。
6. **TED/ECOD**:TED/ECOD 的坐标来自 AFdb v4。只有 AFdb v4 序列与所选结构序列完全一致时
   才纳入 TED/ECOD(当前 4,950 个基因);v4/v6 不一致的基因不做 TED/ECOD,避免坐标错位。
{dy_block}"""


def render_report(
    shared_comparison: pd.DataFrame,
    afdb_only: pd.DataFrame,
    pombase_only: pd.DataFrame,
    counts: dict[str, int],
    config: Config,
    dy_result: Dy47073Result | None = None,
) -> str:
    """Assemble the complete Markdown report."""
    dy_header = ""
    title = "# AFdb v6 与 PomBase 蛋白注释差异报告"
    if dy_result is not None:
        dy_header = f"DY47073 数据:`{config.dy47073_fasta}`(release `{config.dy47073_label}`)\n"
        title = "# AFdb v6、PomBase 与 DY47073 蛋白注释差异报告"
    header = f"""{title}

生成日期:{date.today().isoformat()}
AFdb 数据:`{config.afdb_dir}`
PomBase 数据:`{config.pombase_dir}`(release `{config.pombase_release}`)
{dy_header}
"""
    sections = [
        render_set_overview(counts, config),
        render_length_section(shared_comparison, counts),
        render_sequence_section(shared_comparison, counts),
        render_afdb_only_section(afdb_only, counts),
        render_pombase_only_section(pombase_only, counts),
    ]
    if dy_result is not None:
        sections.append(render_dy47073_section(dy_result, counts, config))
    sections.append(render_recommendations(counts, dy_result))
    return header + "\n".join(sections)


# =============================================================================
# 7. IMPERATIVE SHELL (CLI, Orchestration, Exit Code Propagation)
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--afdb-dir", type=Path, required=True, help="AFdb CIF directory.")
    parser.add_argument("--pae-dir", type=Path, required=True, help="AFdb PAE JSON directory.")
    parser.add_argument(
        "--pombase-dir", type=Path, required=True, help="PomBase release directory."
    )
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument(
        "--pombase-release",
        default="2026-11-01",
        help="PomBase release label written to the report.",
    )
    parser.add_argument(
        "--dy47073-fasta",
        type=Path,
        default=None,
        help="Optional DY47073 protein FASTA for the three-way comparison.",
    )
    parser.add_argument(
        "--dy47073-label",
        default="240717",
        help="DY47073 release label written to the report.",
    )
    parser.add_argument(
        "--uniprot-resolution",
        type=Path,
        default=None,
        help="Optional UniProt accession resolution TSV.",
    )
    parser.add_argument(
        "--uniprot-triangulation",
        type=Path,
        default=None,
        help="Optional AFdb/PomBase/current-UniProt triangulation TSV.",
    )
    return parser.parse_args()


def build_config(args: argparse.Namespace) -> Config:
    """Build the immutable runtime configuration from parsed CLI arguments."""
    return Config(
        afdb_dir=args.afdb_dir.resolve(),
        pae_dir=args.pae_dir.resolve(),
        pombase_dir=args.pombase_dir.resolve(),
        output_dir=args.output_dir.resolve(),
        pombase_release=args.pombase_release,
        uniprot_resolution=(args.uniprot_resolution.resolve() if args.uniprot_resolution else None),
        uniprot_triangulation=(
            args.uniprot_triangulation.resolve() if args.uniprot_triangulation else None
        ),
        dy47073_fasta=(args.dy47073_fasta.resolve() if args.dy47073_fasta else None),
        dy47073_label=args.dy47073_label,
    )


def write_tables(
    config: Config,
    shared_comparison: pd.DataFrame,
    shared_by_sequence: pd.DataFrame,
    merged: pd.DataFrame,
    afdb_entry_table: pd.DataFrame,
    afdb_only: pd.DataFrame,
    pombase_only: pd.DataFrame,
    counts: dict[str, int],
    dy_result: Dy47073Result | None = None,
) -> None:
    """Write all tabular outputs."""
    config.output_dir.mkdir(parents=True, exist_ok=True)

    nonidentical = shared_comparison.loc[shared_comparison[ComparisonCol.CATEGORY] != "identical"]
    nonidentical.to_csv(config.output_dir / "nonidentical_pairs.tsv", sep="\t", index=False)
    shared_comparison.to_csv(config.output_dir / "shared_comparison.tsv", sep="\t", index=False)
    shared_by_sequence.to_csv(config.output_dir / "shared_by_sequence.tsv", sep="\t", index=False)
    merged.to_csv(config.output_dir / "merged_table.tsv", sep="\t", index=False)
    afdb_entry_table.to_csv(config.output_dir / "afdb_entries.tsv", sep="\t", index=False)
    afdb_only.to_csv(config.output_dir / "afdb_only.tsv", sep="\t", index=False)
    pombase_only.to_csv(config.output_dir / "pombase_only.tsv", sep="\t", index=False)
    if dy_result is not None:
        comparison_columns = [*DY47073_COMPARISON_COLUMNS, DyCol.SEQUENCE.value]
        dy_result.comparison.loc[:, comparison_columns].to_csv(
            config.output_dir / "dy47073_comparison.tsv", sep="\t", index=False
        )
        dy_result.only.loc[:, DY47073_ONLY_COLUMNS].to_csv(
            config.output_dir / "dy47073_only.tsv", sep="\t", index=False
        )
        dy_result.afdb_only_matches.loc[:, AFDB_ONLY_DY_COLUMNS].to_csv(
            config.output_dir / "afdb_only_dy47073.tsv", sep="\t", index=False
        )
    pd.DataFrame({"metric": list(counts), "value": list(counts.values())}).to_csv(
        config.output_dir / "summary_counts.tsv", sep="\t", index=False
    )


def write_excel(
    config: Config,
    counts: dict[str, int],
    shared_by_id: pd.DataFrame,
    shared_by_sequence: pd.DataFrame,
    merged: pd.DataFrame,
    afdb_entry_table: pd.DataFrame,
    afdb_only: pd.DataFrame,
    pombase_only: pd.DataFrame,
    dy_result: Dy47073Result | None = None,
) -> None:
    """Write the multi-sheet Excel workbook."""
    summary = pd.DataFrame({"metric": list(counts), "value": list(counts.values())})
    shared = merged.loc[merged["category"] == "shared"]
    workbook_path = config.output_dir / "afdb_vs_pombase.xlsx"

    with pd.ExcelWriter(workbook_path, engine=EXCEL_ENGINE) as writer:
        summary.to_excel(writer, sheet_name="summary", index=False)
        merged.to_excel(writer, sheet_name="merged", index=False)
        afdb_entry_table.to_excel(writer, sheet_name="afdb_entries", index=False)
        shared_by_id.drop(columns=["afdb_sequence", "pombase_sequence"]).to_excel(
            writer, sheet_name="shared_by_id", index=False
        )
        shared_by_sequence.to_excel(writer, sheet_name="shared_by_sequence", index=False)
        if "relationship" in shared_by_id.columns:
            check_columns = [
                column
                for column in shared_by_id.columns
                if column not in ("afdb_sequence", "pombase_sequence")
            ]
            shared_by_id.loc[shared_by_id["relationship"].notna(), check_columns].to_excel(
                writer, sheet_name="shared_uniprot_check", index=False
            )
        for category, sheet_name in (
            ("identical", "shared_identical"),
            ("substitution_only", "shared_substitution"),
            ("indel_only", "shared_indel"),
            ("both", "shared_both"),
        ):
            shared.loc[shared["sequence_category"] == category].to_excel(
                writer, sheet_name=sheet_name, index=False
            )

        afdb_only_no_sequence = afdb_only.drop(columns=["afdb_sequence"])
        afdb_only_no_sequence.to_excel(writer, sheet_name="afdb_only_all", index=False)
        afdb_only_no_sequence.loc[
            afdb_only_no_sequence["match_status"] == "no_pombase_match"
        ].to_excel(writer, sheet_name="afdb_only_unmatched", index=False)
        afdb_only_no_sequence.loc[
            afdb_only_no_sequence["match_status"] == "sequence_identical_to_shared_protein"
        ].to_excel(writer, sheet_name="afdb_only_duplicates", index=False)

        pombase_only.to_excel(writer, sheet_name="pombase_only_all", index=False)
        pombase_only.loc[
            pombase_only["reason"] != "represented_in_afdb_under_other_accession"
        ].to_excel(writer, sheet_name="pombase_only_unrepresented", index=False)

        if dy_result is not None:
            dy_comparison = dy_result.comparison.loc[:, DY47073_COMPARISON_COLUMNS]
            dy_comparison.to_excel(writer, sheet_name="dy47073_comparison", index=False)
            dy_comparison.loc[dy_comparison[DyCol.CATEGORY] != "identical"].to_excel(
                writer, sheet_name="dy47073_differences", index=False
            )
            dy_result.only.loc[:, DY47073_ONLY_COLUMNS].to_excel(
                writer, sheet_name="dy47073_only", index=False
            )
            dy_result.afdb_only_matches.loc[:, AFDB_ONLY_DY_COLUMNS].to_excel(
                writer, sheet_name="afdb_only_dy47073", index=False
            )
            if not dy_result.secondary_matches.empty:
                secondary_columns = [
                    DyCol.GENE,
                    DyCol.RECORD_ID,
                    DyCol.LENGTH,
                    AfdbCol.ACCESSION,
                    AfdbCol.LENGTH,
                    DyCol.VS_POMBASE,
                ]
                dy_result.secondary_matches.loc[:, secondary_columns].to_excel(
                    writer, sheet_name="dy47073_secondary", index=False
                )

    logger.info(f"Excel workbook written to {workbook_path}")


def build_dy47073_result(
    dy47073_fasta: Path,
    gene_catalog: pd.DataFrame,
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    shared_by_id: pd.DataFrame,
    afdb_only: pd.DataFrame,
) -> Dy47073Result:
    """Build the DY47073 comparison tables from validated inputs."""
    proteins = load_dy47073_proteins(dy47073_fasta, gene_catalog)
    afdb_sequences = build_afdb_sequence_table(shared_by_id, afdb_only)
    genes = build_dy47073_gene_table(proteins, pombase_peptides, afdb_sequences)
    comparison = build_dy47073_comparison(pombase_genes, pombase_peptides, afdb_sequences, genes)
    has_triangulation = "relationship" in shared_by_id.columns
    comparison = (
        comparison.merge(
            shared_by_id.loc[:, [GeneCol.SYSTEMATIC_ID, "relationship"]],
            on=GeneCol.SYSTEMATIC_ID,
            how="left",
        )
        if has_triangulation
        else comparison.assign(relationship=np.nan)
    )
    return Dy47073Result(
        comparison=comparison,
        only=build_dy47073_only_table(genes, proteins, gene_catalog, afdb_only),
        afdb_only_matches=annotate_afdb_only_with_dy47073(afdb_only, proteins),
        secondary_matches=build_dy47073_secondary_matches(proteins, afdb_only, comparison),
        record_count=len(proteins),
        internal_stop_count=int(proteins[DyCol.INTERNAL_STOP].sum()),
    )


def main() -> int:
    """Load inputs, run the comparison, and write the report."""
    args = parse_args()
    config = build_config(args)

    try:
        afdb_entries = load_afdb_entries(config.afdb_dir)
        pae_accessions = load_pae_accessions(config.pae_dir)
        gene_catalog = load_pombase_gene_catalog(config.pombase_dir)
        pombase_genes = gene_catalog.loc[
            gene_catalog[GeneCol.GENE_TYPE] == PROTEIN_CODING
        ].reset_index(drop=True)
        pombase_peptides = load_pombase_peptides(config.pombase_dir)

        shared_accessions, afdb_only_accessions, pombase_only_accessions = split_accession_sets(
            afdb_entries, pombase_genes
        )
        shared_genes = pombase_genes.loc[pombase_genes[GeneCol.ACCESSION].isin(shared_accessions)]
        shared_comparison = attach_comparison(
            shared_genes,
            afdb_entries.loc[afdb_entries[AfdbCol.ACCESSION].isin(shared_accessions)],
            pombase_peptides,
            build_aligner(),
        )
        if config.uniprot_triangulation is not None:
            shared_comparison = shared_comparison.merge(
                load_uniprot_triangulation(config.uniprot_triangulation),
                on=[GeneCol.SYSTEMATIC_ID, ComparisonCol.ACCESSION],
                how="left",
            )
        afdb_only = annotate_afdb_only(
            afdb_entries.loc[afdb_entries[AfdbCol.ACCESSION].isin(afdb_only_accessions)],
            pombase_peptides,
            afdb_entries.loc[afdb_entries[AfdbCol.ACCESSION].isin(shared_accessions)],
            shared_genes,
        )
        if config.uniprot_resolution is not None:
            pombase_accession_set = set(pombase_genes[GeneCol.ACCESSION].dropna())
            afdb_only = afdb_only.merge(
                load_uniprot_resolution(config.uniprot_resolution, pombase_accession_set),
                on=AfdbCol.ACCESSION,
                how="left",
            )
        afdb_max_length = int(afdb_entries[AfdbCol.LENGTH].max())
        pombase_only = pd.concat(
            [
                build_pombase_only_missing_accessions(
                    pombase_genes,
                    pombase_peptides,
                    pombase_only_accessions,
                    afdb_only,
                    afdb_max_length,
                ),
                build_pombase_only(pombase_genes, pombase_peptides, afdb_only, afdb_max_length),
            ],
            ignore_index=True,
        )
        counts = build_summary_counts(
            afdb_entries,
            pae_accessions,
            pombase_genes,
            pombase_peptides,
            shared_accessions,
            afdb_only,
            pombase_only_accessions,
            shared_comparison,
            pombase_only,
        )
        shared_by_id = add_shared_match_columns(shared_comparison, afdb_only)
        shared_by_sequence = build_shared_by_sequence(afdb_only, pombase_genes, pombase_peptides)
        merged = build_merged_table(
            shared_by_id,
            shared_by_sequence,
            build_afdb_only_rows(afdb_only),
            build_pombase_only_rows(pombase_only),
        )
        afdb_entry_table = build_afdb_entry_table(
            afdb_entries, pae_accessions, shared_by_id, shared_by_sequence, afdb_only
        )
        counts.update(uniprot_triangulation_counts(shared_comparison))
        counts.update(
            {
                "shared_by_id_genes": len(shared_by_id),
                "shared_by_sequence_genes": len(shared_by_sequence),
                "merged_shared_rows": len(shared_by_id) + len(shared_by_sequence),
                "afdb_only_proteins": int((merged["category"] == "afdb_only").sum()),
                "pombase_only_proteins": int((merged["category"] == "pombase_only").sum()),
                "merged_rows": len(merged),
            }
        )

        dy_result = None
        if config.dy47073_fasta is not None:
            dy_result = build_dy47073_result(
                config.dy47073_fasta,
                gene_catalog,
                pombase_genes,
                pombase_peptides,
                shared_by_id,
                afdb_only,
            )
            merged = attach_dy47073_to_merged(merged, dy_result)
            counts.update(dy47073_counts(dy_result))

        write_tables(
            config,
            shared_comparison,
            shared_by_sequence,
            merged,
            afdb_entry_table,
            afdb_only,
            pombase_only,
            counts,
            dy_result,
        )
        write_excel(
            config,
            counts,
            shared_by_id,
            shared_by_sequence,
            merged,
            afdb_entry_table,
            afdb_only,
            pombase_only,
            dy_result,
        )
        report = render_report(
            shared_comparison, afdb_only, pombase_only, counts, config, dy_result
        )
        (config.output_dir / "report.md").write_text(report, encoding="utf-8")
    except (OSError, ValueError, KeyError) as error:
        logger.exception(f"Fatal execution error: {error}")
        return 1
    else:
        logger.info(f"Report written to {config.output_dir / 'report.md'}")
        return 0


if __name__ == "__main__":
    sys.exit(main())
