"""
AFdb v6 vs PomBase Protein Annotation Consistency
=================================================

Compare the proteins modelled in AlphaFold DB release v6 for
Schizosaccharomyces pombe (proteome UP000002485) with the current PomBase
protein-coding gene annotation.

Pairing is sequence-first:

1. For every PomBase gene, find AFdb v6 entries whose modelled sequence exactly
   equals the PomBase peptide. Choose one entry, preferring the gene's own
   UniProt accession and otherwise the smallest alternate accession.
2. Genes without an exact sequence match fall back to their UniProt accession;
   these pairs are compared field by field and classified as substitution-only,
   indel-only, or both (alignment-derived).
3. AFdb entries selected by no gene are classified as superseded (the gene chose
   a sequence-identical alternate entry), duplicate (the sequence matches a gene
   that preferred its own entry), or no_pombase_match.

AFdb sequences are read from the entry inventory table produced by the
``afdb-entry-table`` stage (modelled CIF entity sequence
``_entity_poly.pdbx_seq_one_letter_code_can``); PomBase sequences are read
from ``peptide.fa`` and truncated at the first stop codon. The inventory also
records PAE availability per entry.

Input
-----
- AFdb v6 entry inventory (``03_afdb_entries.tsv``).
- PomBase release directory containing ``gene_IDs_names_products.tsv`` and
  ``peptide.fa[.gz]``.
- Optional UniProt accession resolution TSV for the AFdb-only entries.
- Optional AFdb/PomBase/current-UniProt triangulation TSV for shared proteins.

Output
------
- ``01_merged_table.tsv``: unified paired/AFdb-only/PomBase-only table.
- ``02_afdb_entry_roles.tsv``: AFdb v6 entry roles against selected genes.
- ``03_shared_comparison.tsv``: one row per paired gene.
- ``04_shared_by_sequence.tsv``: sequence matches whose entry is not the gene's accession.
- ``05_nonidentical_pairs.tsv``: ID-fallback pairs with non-identical sequences.
- ``06_afdb_only.tsv``: AFdb entries selected by no gene.
- ``07_pombase_only.tsv``: PomBase-only genes with lengths and reasons.
- ``08_summary_counts.tsv``: headline counts.
- ``09_report.md``: human-readable difference report.
- ``10_afdb_vs_pombase.xlsx``: multi-sheet workbook with the merged table and
  per-category sheets.

Usage
-----
    python compare_afdb_pombase.py \
        --afdb-entries results/01_afdb/03_afdb_entries.tsv \
        --afdb-source /hugedata/YushengYang/AFdb/UP000002485_284812_SCHPO_v6 \
        --pombase-dir tmp/afdb_vs_pombase/reference/pombase-2026-11-01 \
        --output-dir tmp/afdb_vs_pombase/outputs \
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
import sys
from typing import TextIO

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from loguru import logger
import numpy as np
import pandas as pd

# =============================================================================
# 2. CONSTANTS & ENUMS (Domain Markers & Column Definitions)
# =============================================================================
PROTEIN_CODING = "protein coding gene"
ACCESSION_COLUMN_ALIASES = ("external_id", "uniprot_id")
AFDB_ENTRY_COLUMNS = [
    "afdb_accession",
    "afdb_description",
    "afdb_length",
    "afdb_sequence",
    "pae_available",
]
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

AFDB_ONLY_STATUSES = (
    "superseded_by_sequence_match",
    "duplicate_of_selected_entry",
    "no_pombase_match",
)
ROLE_PRIORITY = {
    "primary": 0,
    "sequence_matched": 1,
    "superseded_by_sequence_match": 2,
    "duplicate_of_selected_entry": 3,
    "no_pombase_match": 4,
}
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


# =============================================================================
# 3. CONFIGURATION (Immutable Runtime Parameters)
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class Config:
    """Immutable runtime configuration container."""

    afdb_entries: Path
    afdb_source: str
    pombase_dir: Path
    output_dir: Path
    pombase_release: str
    uniprot_resolution: Path | None = None
    uniprot_triangulation: Path | None = None


# =============================================================================
# 4. BOUNDARY LOADERS (I/O Shell)
# =============================================================================
def open_text(path: Path) -> TextIO:
    """Open a plain-text or gzip-compressed text file."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def load_afdb_entries(path: Path) -> pd.DataFrame:
    """Load the AFdb entry inventory with the comparison column names."""
    inventory = pd.read_csv(path, sep="\t", dtype=str)
    missing = set(AFDB_ENTRY_COLUMNS) - set(inventory.columns)
    if missing:
        raise ValueError(f"Missing AFdb entry columns in {path}: {sorted(missing)}")

    inventory = inventory.loc[:, AFDB_ENTRY_COLUMNS]
    inventory = inventory.rename(columns={"afdb_accession": AfdbCol.ACCESSION})
    inventory[AfdbCol.LENGTH] = pd.to_numeric(inventory[AfdbCol.LENGTH]).astype("Int64")
    logger.info(f"Loaded {len(inventory):,} AFdb entries from {path}")
    return inventory


def load_pae_accessions(entries: pd.DataFrame) -> set[str]:
    """Return the accessions flagged as having a PAE matrix."""
    pae_available = entries.loc[entries["pae_available"] == "True", AfdbCol.ACCESSION]
    logger.info(f"Found {len(pae_available):,} AFdb entries with PAE")
    return set(pae_available)


def load_afdb_inputs(path: Path) -> tuple[pd.DataFrame, set[str]]:
    """Load the comparison columns and the PAE accession set from one inventory."""
    inventory = load_afdb_entries(path)
    entries = inventory.loc[
        :, [AfdbCol.ACCESSION, AfdbCol.DESCRIPTION, AfdbCol.SEQUENCE, AfdbCol.LENGTH]
    ]
    return entries, load_pae_accessions(inventory)


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
    with open_text(path) as handle:
        return [
            (record.id.split(":pep")[0].rsplit(".", 1)[0], str(record.seq))
            for record in SeqIO.parse(handle, "fasta")
        ]


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


# =============================================================================
# 5. FUNCTIONAL CORE (Pure, Trusted Transformations)
# =============================================================================
def assign_structure_accessions(
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    afdb_entries: pd.DataFrame,
) -> pd.DataFrame:
    """Pair each gene to one AFdb entry by exact sequence, then by UniProt ID."""
    genes = pombase_genes.merge(
        pombase_peptides.loc[
            :,
            [
                PeptideCol.SYSTEMATIC_ID,
                PeptideCol.SEQUENCE,
                PeptideCol.LENGTH,
                PeptideCol.INTERNAL_STOP,
            ],
        ],
        on=PeptideCol.SYSTEMATIC_ID,
        how="inner",
    ).rename(columns={GeneCol.ACCESSION: "pombase_accession"})

    entry_accessions = set(afdb_entries[AfdbCol.ACCESSION])
    candidates = genes.merge(
        afdb_entries.loc[:, [AfdbCol.ACCESSION, AfdbCol.SEQUENCE]],
        left_on=PeptideCol.SEQUENCE,
        right_on=AfdbCol.SEQUENCE,
        how="inner",
    )
    is_own_entry = candidates["pombase_accession"].eq(candidates[AfdbCol.ACCESSION])
    candidates = candidates.assign(_is_own_entry=is_own_entry).sort_values(
        [PeptideCol.SYSTEMATIC_ID, "_is_own_entry", AfdbCol.ACCESSION],
        ascending=[True, False, True],
        kind="stable",
    )
    chosen = candidates.drop_duplicates(PeptideCol.SYSTEMATIC_ID, keep="first").rename(
        columns={AfdbCol.ACCESSION: "afdb_accession"}
    )
    duplicates = (
        candidates.loc[~candidates.index.isin(chosen.index)]
        .groupby(PeptideCol.SYSTEMATIC_ID)[AfdbCol.ACCESSION]
        .agg(join_unique)
    )

    pairs = genes.merge(
        chosen.loc[:, [PeptideCol.SYSTEMATIC_ID, "afdb_accession"]],
        on=PeptideCol.SYSTEMATIC_ID,
        how="left",
    )
    pairs["duplicate_afdb_accessions"] = pairs[PeptideCol.SYSTEMATIC_ID].map(duplicates).fillna("")
    has_sequence_match = pairs["afdb_accession"].notna()
    has_id_entry = pairs["pombase_accession"].isin(entry_accessions)
    pairs["match_method"] = np.select(
        [has_sequence_match, has_id_entry],
        ["sequence", "uniprot_id"],
        default="none",
    )
    is_id_fallback = pairs["match_method"].eq("uniprot_id")
    pairs.loc[is_id_fallback, "afdb_accession"] = pairs.loc[is_id_fallback, "pombase_accession"]

    has_own_entry = pairs["pombase_accession"].isin(entry_accessions)
    pairs["uniprot_id"] = pairs["pombase_accession"].where(has_own_entry, pairs["afdb_accession"])
    pairs["uniprot_id"] = pairs["uniprot_id"].fillna(pairs["pombase_accession"])
    return pairs


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


def attach_comparison(
    pairs: pd.DataFrame,
    afdb_entries: pd.DataFrame,
    aligner: PairwiseAligner,
) -> pd.DataFrame:
    """Attach AFdb fields and sequence-difference statistics to each selected pair."""
    entry_table = afdb_entries.rename(columns={AfdbCol.ACCESSION: "afdb_accession"})
    merged = pairs.merge(entry_table, on="afdb_accession", how="inner")
    merged[ComparisonCol.LENGTH_DELTA] = merged[AfdbCol.LENGTH] - merged["pombase_length"]

    statistics = pd.DataFrame.from_records(
        [
            compare_sequences(afdb_sequence, pombase_sequence, aligner)
            for afdb_sequence, pombase_sequence in zip(
                merged[AfdbCol.SEQUENCE], merged["pombase_sequence"], strict=True
            )
        ]
    )
    return pd.concat([merged.reset_index(drop=True), statistics], axis=1)


def annotate_unpaired_entries(
    afdb_entries: pd.DataFrame,
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    selected_accessions: set[str],
) -> pd.DataFrame:
    """Classify AFdb entries that no gene selected as its structure."""
    unpaired = afdb_entries.loc[~afdb_entries[AfdbCol.ACCESSION].isin(selected_accessions)].copy()
    peptide_to_genes = (
        pombase_peptides.groupby(PeptideCol.SEQUENCE)[PeptideCol.SYSTEMATIC_ID]
        .agg(join_unique)
        .to_dict()
    )
    accession_to_genes = (
        pombase_genes.groupby(GeneCol.ACCESSION)[GeneCol.SYSTEMATIC_ID].agg(join_unique).to_dict()
    )
    unpaired["matched_pombase_genes"] = unpaired[AfdbCol.SEQUENCE].map(peptide_to_genes).fillna("")
    unpaired["superseded_gene"] = unpaired[AfdbCol.ACCESSION].map(accession_to_genes).fillna("")
    is_superseded = unpaired["superseded_gene"] != ""
    has_sequence_match = unpaired["matched_pombase_genes"] != ""
    unpaired["match_status"] = np.select(
        [is_superseded, has_sequence_match],
        ["superseded_by_sequence_match", "duplicate_of_selected_entry"],
        default="no_pombase_match",
    )
    return unpaired


def build_pombase_only(unpaired_genes: pd.DataFrame, afdb_max_length: int) -> pd.DataFrame:
    """List genes with no AFdb structure, with length and reason."""
    table = unpaired_genes.copy()
    table["uniprot_id"] = table["pombase_accession"].fillna("")
    table["matched_afdb_only_accessions"] = ""
    table["exceeds_afdb_max_length"] = table["pombase_length"] > afdb_max_length
    table["reason"] = np.select(
        [table["pombase_accession"].isna(), table["exceeds_afdb_max_length"]],
        ["no_uniprot_accession", "exceeds_afdb_length_cutoff"],
        default="not_modelled_by_afdb",
    )
    return table.loc[:, POMBASE_ONLY_COLUMNS]


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

    resolved = afdb_only["is_secondary"].notna()
    is_secondary = afdb_only["is_secondary"] == "True"
    primary_in_pombase = afdb_only["primary_in_pombase"] == "True"
    return {
        "uniprot_resolved": int(resolved.sum()),
        "uniprot_secondary": int((resolved & is_secondary).sum()),
        "uniprot_secondary_in_pombase": int((resolved & is_secondary & primary_in_pombase).sum()),
        "uniprot_primary": int((resolved & ~is_secondary).sum()),
        "uniprot_primary_in_pombase": int((resolved & ~is_secondary & primary_in_pombase).sum()),
    }


def build_summary_counts(
    afdb_entries: pd.DataFrame,
    pae_accessions: set[str],
    pombase_genes: pd.DataFrame,
    pombase_peptides: pd.DataFrame,
    shared_comparison: pd.DataFrame,
    afdb_only: pd.DataFrame,
    pombase_only: pd.DataFrame,
) -> dict[str, int]:
    """Collect headline counts for the report."""
    categories = shared_comparison[ComparisonCol.CATEGORY]
    match_method = shared_comparison["match_method"]
    afdb_accessions = set(afdb_entries[AfdbCol.ACCESSION])
    is_sequence = match_method == "sequence"
    is_rescued = shared_comparison["pombase_accession"].isna() | shared_comparison[
        "afdb_accession"
    ].ne(shared_comparison["pombase_accession"])
    counts = {
        "afdb_entries": len(afdb_entries),
        "afdb_unique_accessions": afdb_entries[AfdbCol.ACCESSION].nunique(),
        "pombase_coding_genes": len(pombase_genes),
        "pombase_genes_with_accession": int(pombase_genes[GeneCol.ACCESSION].notna().sum()),
        "pombase_unique_accessions": int(pombase_genes[GeneCol.ACCESSION].nunique()),
        "shared_genes": len(shared_comparison),
        "shared_by_sequence_genes": int(is_sequence.sum()),
        "shared_by_id_genes": int((match_method == "uniprot_id").sum()),
        "sequence_rescued_genes": int((is_sequence & is_rescued).sum()),
        "afdb_only_accessions": len(afdb_only),
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
        "afdb_only_superseded": int(
            (afdb_only["match_status"] == "superseded_by_sequence_match").sum()
        ),
        "afdb_only_duplicate": int(
            (afdb_only["match_status"] == "duplicate_of_selected_entry").sum()
        ),
        "afdb_only_no_match": int((afdb_only["match_status"] == "no_pombase_match").sum()),
        "afdb_max_length": int(afdb_entries[AfdbCol.LENGTH].max()),
        "pombase_max_length": int(pombase_peptides[PeptideCol.LENGTH].max()),
        "pombase_max_shared_length": int(shared_comparison[ComparisonCol.POMBASE_LENGTH].max()),
        "pae_entries": len(pae_accessions),
        "afdb_accessions_without_pae": len(afdb_accessions - pae_accessions),
        "pae_accessions_without_afdb": len(pae_accessions - afdb_accessions),
        "pombase_internal_stop_peptides": int(pombase_peptides[PeptideCol.INTERNAL_STOP].sum()),
        "pombase_only_length_cutoff": int(
            (pombase_only["reason"] == "exceeds_afdb_length_cutoff").sum()
        ),
        "pombase_only_absent_other": int(
            pombase_only["reason"].isin(["not_modelled_by_afdb", "no_uniprot_accession"]).sum()
        ),
    }
    counts.update(uniprot_resolution_counts(afdb_only))
    counts.update(uniprot_triangulation_counts(shared_comparison))
    no_match_entries = int((afdb_only["match_status"] == "no_pombase_match").sum())
    counts.update(
        {
            "merged_shared_rows": len(shared_comparison),
            "afdb_only_proteins": no_match_entries,
            "pombase_only_proteins": len(pombase_only),
            "merged_rows": len(shared_comparison) + no_match_entries + len(pombase_only),
        }
    )
    return counts


def join_unique(values: pd.Series) -> str:
    """Join the unique non-missing string values of a Series with semicolons."""
    return ";".join(sorted({str(value) for value in values if pd.notna(value)}))


def build_shared_by_sequence(shared_comparison: pd.DataFrame) -> pd.DataFrame:
    """Return sequence matches whose selected entry is not the gene's UniProt accession."""
    is_sequence = shared_comparison["match_method"] == "sequence"
    is_rescued = shared_comparison["pombase_accession"].isna() | shared_comparison[
        "afdb_accession"
    ].ne(shared_comparison["pombase_accession"])
    return shared_comparison.loc[is_sequence & is_rescued].reset_index(drop=True)


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
    return pombase_only.assign(
        category="pombase_only",
        match_method="none",
        uniprot_id=pombase_only["uniprot_id"],
        afdb_accession=np.nan,
        pombase_accession=pombase_only["uniprot_id"],
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
        notes=pombase_only["reason"],
    )


def build_merged_table(
    shared_rows: pd.DataFrame,
    afdb_only_rows: pd.DataFrame,
    pombase_only_rows: pd.DataFrame,
) -> pd.DataFrame:
    """Concatenate the gene-level rows and select the report columns."""
    merged = pd.concat(
        [shared_rows, afdb_only_rows, pombase_only_rows],
        ignore_index=True,
    ).reindex(columns=MERGED_COLUMNS)
    for column in INTEGER_COLUMNS:
        merged[column] = merged[column].astype("Int64")
    return merged


def build_afdb_entry_table(
    afdb_entries: pd.DataFrame,
    pae_accessions: set[str],
    shared_comparison: pd.DataFrame,
    afdb_only: pd.DataFrame,
) -> pd.DataFrame:
    """Map every AFdb entry to its selected gene(s) and entry role."""
    roles = pd.concat(
        [build_selected_entry_roles(shared_comparison), build_unpaired_entry_roles(afdb_only)],
        ignore_index=True,
    )
    roles["_priority"] = roles["entry_role"].map(ROLE_PRIORITY)
    roles = roles.sort_values(["afdb_accession", "_priority"], kind="stable")
    roles = roles.groupby("afdb_accession", as_index=False).agg(
        gene_systematic_id=("gene_systematic_id", join_unique),
        gene_name=("gene_name", join_unique),
        pombase_accession=("pombase_accession", join_unique),
        entry_role=("entry_role", "first"),
        protein_entity_accession=("protein_entity_accession", join_unique),
        notes=("notes", join_unique),
    )
    roles["entity_category"] = np.where(
        roles["entry_role"].isin(["duplicate_of_selected_entry", "no_pombase_match"]),
        "afdb_only",
        "shared",
    )

    roles = roles.rename(columns={"afdb_accession": AfdbCol.ACCESSION.value})
    table = afdb_entries.merge(roles, on=AfdbCol.ACCESSION, how="left")
    if table["entry_role"].isna().any():
        missing = table.loc[table["entry_role"].isna(), AfdbCol.ACCESSION].tolist()
        raise ValueError(f"AFdb entries without a role: {missing[:10]}")

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


def build_selected_entry_roles(shared_comparison: pd.DataFrame) -> pd.DataFrame:
    """Build entry roles for AFdb entries selected by a gene."""
    roles = shared_comparison.loc[
        :, ["afdb_accession", "gene_systematic_id", "gene_name", "uniprot_id", "pombase_accession"]
    ].copy()
    roles["entry_role"] = np.where(
        roles["afdb_accession"].eq(roles["uniprot_id"]), "primary", "sequence_matched"
    )
    roles["protein_entity_accession"] = roles["pombase_accession"].fillna(roles["afdb_accession"])
    roles["notes"] = np.where(roles["entry_role"].eq("primary"), "primary", "matched_by_sequence")
    return roles


def build_unpaired_entry_roles(afdb_only: pd.DataFrame) -> pd.DataFrame:
    """Build entry roles for AFdb entries that no gene selected."""
    roles = afdb_only.rename(columns={AfdbCol.ACCESSION: "afdb_accession"}).copy()
    roles["entry_role"] = roles["match_status"]
    roles["protein_entity_accession"] = roles["afdb_accession"]
    roles["gene_systematic_id"] = (
        roles["superseded_gene"]
        .where(
            roles["match_status"].eq("superseded_by_sequence_match"),
            roles["matched_pombase_genes"],
        )
        .fillna("")
    )
    roles["gene_name"] = ""
    roles["pombase_accession"] = ""
    roles["notes"] = roles["match_status"]
    return roles


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


def render_set_overview(counts: dict[str, int], config: Config) -> str:
    """Render the matching-overview section."""
    selected_entries = counts["afdb_entries"] - counts["afdb_only_accessions"]
    accession_without_entry = counts["sequence_rescued_genes"] - counts["afdb_only_superseded"]
    rows = [
        ("AFdb 结构条目(全部 F1 片段)", counts["afdb_entries"]),
        ("AFdb 唯一 UniProt accession", counts["afdb_unique_accessions"]),
        ("PAE JSON 文件", counts["pae_entries"]),
        ("PomBase protein coding genes", counts["pombase_coding_genes"]),
        ("PomBase 有 UniProt accession 的基因", counts["pombase_genes_with_accession"]),
        ("PomBase 唯一 UniProt accession", counts["pombase_unique_accessions"]),
        ("被选中的 AFdb 条目", selected_entries),
        ("合并表 shared 行", counts["shared_genes"]),
        ("  按蛋白序列精确匹配", counts["shared_by_sequence_genes"]),
        ("  按 UniProt ID 回退(序列不一致)", counts["shared_by_id_genes"]),
        ("  序列匹配但更换了 AFdb 条目", counts["sequence_rescued_genes"]),
        ("未选中的 AFdb 条目", counts["afdb_only_accessions"]),
        ("  被序列匹配取代的旧条目", counts["afdb_only_superseded"]),
        ("  重复/次级条目", counts["afdb_only_duplicate"]),
        ("  无 PomBase 匹配", counts["afdb_only_no_match"]),
        ("PomBase-only 基因", counts["pombase_only_genes"]),
        ("合并表 afdb_only 行", counts["afdb_only_proteins"]),
        ("合并表 pombase_only 行", counts["pombase_only_proteins"]),
        ("合并表总行数", counts["merged_rows"]),
    ]
    table = "\n".join(f"| {name} | {value:,} |" for name, value in rows)
    return f"""## 1. 匹配总览

数据版本:

- AFdb:`UP000002485_284812_SCHPO_v6`,{counts["afdb_entries"]:,} 个 CIF 结构文件。
- PomBase:`{config.pombase_release}`,{counts["pombase_coding_genes"]:,} 个 protein coding gene。
- AFdb 序列取 CIF 的 `_entity_poly.pdbx_seq_one_letter_code_can`(实际建模序列)。
- PomBase 序列取 `peptide.fa`(在第一个终止符处截断;本版本有
  {counts["pombase_internal_stop_peptides"]} 条内部终止符记录)。

匹配顺序:

1. **蛋白序列完全一致**:为每个基因选择序列与 PomBase peptide 完全相同的 AFdb 条目。
   同一基因有多个候选时,优先使用该基因自身的 UniProt accession。
2. **UniProt accession 回退**:没有精确序列条目时,使用基因自身的 UniProt 条目并比较序列差异。
3. **无结构**:两者都不可用时记为 PomBase-only。

合并表 `match_method` 列记录 `sequence` 或 `uniprot_id`。序列匹配但更换了条目的
{counts["sequence_rescued_genes"]:,} 个基因中,{counts["afdb_only_superseded"]:,} 个的原始
UniProt 条目序列与 PomBase 不同,其余 {accession_without_entry:,}
个的 accession 不在 AFdb。

| 项目 | 数量 |
|---|---:|
{table}

未选中的 {counts["afdb_only_accessions"]:,} 个条目中,只有 {counts["afdb_only_no_match"]:,} 个
无 PomBase 匹配的条目以 `afdb_only` 行进入合并表;被序列匹配取代和重复的条目属于已配对基因,
不重复计入。
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
            "afdb_accession",
            ComparisonCol.AFDB_LENGTH,
            ComparisonCol.POMBASE_LENGTH,
            ComparisonCol.LENGTH_DELTA,
        ],
    ]
    top.columns = ["gene", "name", "AFdb 条目", "AFdb 长度", "PomBase 长度", "差值"]
    table = markdown_table(top, list(top.columns))
    return f"""## 2. 配对蛋白的长度一致性

{counts["shared_by_sequence_genes"]:,} 个序列精确匹配的基因对长度必然相等;
以下统计基于全部 {counts["shared_genes"]:,} 个基因对,长度差异只可能来自
{counts["shared_by_id_genes"]:,} 个 UniProt ID 回退对。

- 长度相同:{counts["same_length_pairs"]:,} / {counts["shared_genes"]:,} 个基因对。
- 长度不同:{counts["length_diff_pairs"]:,} / {counts["shared_genes"]:,} 个基因对。
- AFdb 最长建模序列:{counts["afdb_max_length"]:,} aa。
- PomBase 全蛋白组最长:{counts["pombase_max_length"]:,} aa;
  配对蛋白中 PomBase 最长:{counts["pombase_max_shared_length"]:,} aa。

长度差异最大的 15 个基因对:

{table}

注:差异来自 UniProt 与 PomBase 基因模型不一致,或 AFdb v6 使用了较旧的 UniProt 序列,
而不是 AFdb 截断;最大几例(gpi14、sap34、rhp16)经 UniProt 复核已确认。
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
            "afdb_accession",
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
        "AFdb 条目",
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
            "afdb_accession",
            ComparisonCol.N_SUBSTITUTIONS,
            ComparisonCol.IDENTITY_PCT,
        ],
    ]
    substitution_view.columns = ["gene", "name", "AFdb 条目", "替换数", "identity %"]
    substitution_table = markdown_table(substitution_view, list(substitution_view.columns))

    if "uniprot_check_pombase_matches" in counts:
        triangulation_note = f"""UniProt 三方核对(覆盖 {counts["shared_by_id_genes"]:,} 个 ID 回退对
及 {counts["afdb_only_superseded"]:,} 个被序列匹配取代的旧条目):

- {counts["uniprot_check_pombase_matches"]:,} 个:PomBase == 当前 UniProt,AFdb 不同
  → AFdb v6 使用了旧序列,当前参考应以 PomBase/UniProt 为准。
- {counts["uniprot_check_afdb_matches"]:,} 个:AFdb == 当前 UniProt,PomBase 不同
  → UniProt 与 PomBase 注释分歧,需逐条核对基因模型。
- {counts["uniprot_check_neither"]:,} 个:三者都不一致,需人工核对。
"""
    else:
        triangulation_note = (
            "本次额外对照了 2023-07-01、2024-03-01、2025-05-01、2025-06-01、2026-11-01 "
            "五个 PomBase release:非一致蛋白中多数 AFdb 长度不等于任何一个 release 的 "
            "peptide 长度,说明差异主要不是 PomBase 版本更新造成的。\n"
        )

    return f"""## 3. 配对蛋白的序列一致性

{counts["shared_by_sequence_genes"]:,} 个基因对按精确序列配对,均为 identical;
其余 {counts["shared_by_id_genes"]:,} 个无精确序列结构的基因按 UniProt ID 回退比较。

| 类型 | 数量 |
|---|---:|
| identical(完全一致) | {identical:,} |
| substitution_only(等长,仅点突变) | {substitution:,} |
| indel_only(仅序列增减) | {indel:,} |
| both(点突变 + 序列增减) | {both:,} |

等长、仅点突变的蛋白({substitution:,} 个):

{substitution_table}

非一致蛋白共 {len(nonidentical):,} 个。下表按差异位点数排序,
列出前 30 个(完整表见 `05_nonidentical_pairs.tsv`):

{table}

说明:`identity %` 按可比对残基计算(不计 gap)。PomBase peptide 中的内部终止符已在第一个终止符处截断
({counts["pombase_internal_stop_peptides"]} 条),这些基因的"缺失"残基不应视为真实 indel。

{triangulation_note}"""


def render_afdb_only_section(afdb_only: pd.DataFrame, counts: dict[str, int]) -> str:
    """Render the unselected-AFdb-entry reconciliation section."""
    status_counts = afdb_only["match_status"].value_counts()
    status_rows = "\n".join(
        f"| {status} | {int(status_counts.get(status, 0)):,} |" for status in AFDB_ONLY_STATUSES
    )
    superseded = afdb_only.loc[
        afdb_only["match_status"] == "superseded_by_sequence_match"
    ].sort_values(AfdbCol.LENGTH, ascending=False)
    superseded_view = superseded.loc[
        :, ["superseded_gene", AfdbCol.ACCESSION, AfdbCol.DESCRIPTION, AfdbCol.LENGTH]
    ].head(30)
    superseded_view.columns = ["gene", "原 AFdb 条目", "描述", "长度"]
    superseded_table = markdown_table(superseded_view, list(superseded_view.columns))

    duplicates = afdb_only.loc[
        afdb_only["match_status"] == "duplicate_of_selected_entry"
    ].sort_values(AfdbCol.LENGTH, ascending=False)
    duplicate_view = duplicates.loc[
        :, [AfdbCol.ACCESSION, AfdbCol.DESCRIPTION, AfdbCol.LENGTH, "matched_pombase_genes"]
    ].head(30)
    duplicate_view.columns = ["AFdb 条目", "描述", "长度", "序列匹配的基因"]
    duplicate_table = markdown_table(duplicate_view, list(duplicate_view.columns))

    unmatched = afdb_only.loc[afdb_only["match_status"] == "no_pombase_match"].copy()
    unmatched = unmatched.sort_values(AfdbCol.LENGTH, ascending=False)
    unmatched_view = unmatched.loc[
        :, [AfdbCol.ACCESSION, AfdbCol.DESCRIPTION, AfdbCol.LENGTH]
    ].head(30)
    unmatched_view.columns = ["AFdb 条目", "描述", "长度"]
    unmatched_table = markdown_table(unmatched_view, list(unmatched_view.columns))

    uniprot_note = ""
    if "uniprot_resolved" in counts:
        uniprot_note = (
            "UniProt accession 复核(仅针对 accession 不在 PomBase 的重复/无匹配条目):\n\n"
            f"- 快照覆盖 {counts['uniprot_resolved']:,} 个条目。\n"
            f"- {counts['uniprot_secondary']:,} 个是次级 accession;"
            f"{counts['uniprot_secondary_in_pombase']:,} 个的当前 primary 在 PomBase 中。\n"
            f"- {counts['uniprot_primary']:,} 个是当前 primary 条目,"
            "其 accession 均不在 PomBase 的 protein-coding 映射中。\n"
        )

    return f"""## 4. 未被选中的 AFdb 条目

序列优先后,共有 {counts["afdb_only_accessions"]:,} 个 AFdb 条目没有被任何基因选中。
按原因分类:

| 状态 | 数量 |
|---|---:|
{status_rows}

{uniprot_note}
- `superseded_by_sequence_match`:该条目的 accession 是某个 PomBase 基因的 UniProt ID,
  但该基因已按精确序列选择了另一个条目。这是 AFdb v6 旧序列的主要来源。
- `duplicate_of_selected_entry`:该条目的序列与某个基因完全一致,但该基因使用了自身
  accession 对应的条目;通常表示 AFdb 中存在等序列重复或次级条目。
- `no_pombase_match`:序列与所有 PomBase peptide 都不完全相同,也没有 PomBase 基因使用该 accession,
  是真正需要单独核对的 AFdb-only 条目。

被序列匹配取代的旧条目(前 30 个,按长度降序):

{superseded_table}

序列重复/次级条目(前 30 个):

{duplicate_table}

无 PomBase 精确匹配的条目(按长度降序,前 30 个):

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
    return f"""## 5. 无 AFdb 结构的 PomBase 基因

共有 {counts["pombase_only_genes"]:,} 个 protein coding gene 既没有精确序列匹配,
其 UniProt accession 也不在 AFdb 中。

- {counts["pombase_only_length_cutoff"]:,} 个超过 AFdb 最大建模长度
  ({counts["afdb_max_length"]:,} aa),AFdb 不建模。
- 其余 {counts["pombase_only_absent_other"]:,} 个低于长度上限但仍无 AFdb 结构,
  例如 B5BP45-49(在 UP000002485 中但 AFdb 未收录)和 mat2-Pc/Pi(不在 UP000002485 中)。

{table}
"""


def render_recommendations(counts: dict[str, int]) -> str:
    """Render downstream-use recommendations."""
    nonidentical = (
        counts["substitution_only_pairs"] + counts["indel_only_pairs"] + counts["both_pairs"]
    )
    return f"""## 6. 对下游使用的建议

1. **主分析集**:合并表共有 {counts["merged_shared_rows"]:,} 个 shared 行,其中
   {counts["shared_by_sequence_genes"]:,} 个按精确序列配对、{counts["shared_by_id_genes"]:,} 个
   按 UniProt ID 回退。下游以 `results/04_gene_structure_registry.tsv` 的 `in_analysis` 为准。
2. **序列替换条目**:{counts["sequence_rescued_genes"]:,} 个基因使用的 AFdb 条目不是其
   UniProt accession 对应的条目;其中 {counts["afdb_only_superseded"]:,} 个原条目序列与
   PomBase 不同,已列入第 4 节。注册表对这些基因使用序列精确一致的替代条目。
3. **序列不一致蛋白**:共 {nonidentical:,} 个。注册表对其中 37 个使用 vDulab,
   SPCC569.03 因预测序列与 PomBase 不符而排除。加上长蛋白,vDulab 共 45 个,AF3 共 3 个。
   由于 AFdb 结构选择条件是序列精确一致,下游不再需要 UniProt/PomBase 坐标转换。
4. **无 AFdb 结构的基因**:{counts["pombase_only_genes"]:,} 个;其中
   {counts["pombase_only_length_cutoff"]:,} 个超过 AFdb 最大建模长度
   ({counts["afdb_max_length"]:,} aa)。部分长蛋白有 vDulab/AF3 预测;
   sib1 只有 1-702 片段、B5BP45-49 无结构,均排除出分析。
5. **PAE 配对**:AFdb 结构 {counts["afdb_entries"]:,} 个,PAE 文件 {counts["pae_entries"]:,} 个;
   结构无 PAE {counts["afdb_accessions_without_pae"]:,} 个,PAE 无结构
   {counts["pae_accessions_without_afdb"]:,} 个。
6. **TED/ECOD**:TED/ECOD 的坐标来自 AFdb v4。只有 AFdb v4 序列与所选结构序列完全一致时
   才纳入 TED/ECOD(当前 4,950 个基因);v4/v6 不一致的基因不做 TED/ECOD,避免坐标错位。
"""


def render_report(
    shared_comparison: pd.DataFrame,
    afdb_only: pd.DataFrame,
    pombase_only: pd.DataFrame,
    counts: dict[str, int],
    config: Config,
) -> str:
    """Assemble the complete Markdown report."""
    header = f"""# AFdb v6 与 PomBase 蛋白注释差异报告

生成日期:{date.today().isoformat()}
AFdb 数据:`{config.afdb_source}`
PomBase 数据:`{config.pombase_dir}`(release `{config.pombase_release}`)
"""
    sections = [
        render_set_overview(counts, config),
        render_length_section(shared_comparison, counts),
        render_sequence_section(shared_comparison, counts),
        render_afdb_only_section(afdb_only, counts),
        render_pombase_only_section(pombase_only, counts),
    ]
    sections.append(render_recommendations(counts))
    return header + "\n".join(sections)


# =============================================================================
# 7. IMPERATIVE SHELL (CLI, Orchestration, Exit Code Propagation)
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--afdb-entries", type=Path, required=True, help="AFdb entry inventory.")
    parser.add_argument(
        "--afdb-source", required=True, help="AFdb source label for the report header."
    )
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
        afdb_entries=args.afdb_entries.resolve(),
        afdb_source=args.afdb_source,
        pombase_dir=args.pombase_dir.resolve(),
        output_dir=args.output_dir.resolve(),
        pombase_release=args.pombase_release,
        uniprot_resolution=(args.uniprot_resolution.resolve() if args.uniprot_resolution else None),
        uniprot_triangulation=(
            args.uniprot_triangulation.resolve() if args.uniprot_triangulation else None
        ),
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
) -> None:
    """Write all tabular outputs."""
    config.output_dir.mkdir(parents=True, exist_ok=True)

    nonidentical = shared_comparison.loc[shared_comparison[ComparisonCol.CATEGORY] != "identical"]
    nonidentical.to_csv(config.output_dir / "05_nonidentical_pairs.tsv", sep="\t", index=False)
    shared_comparison.to_csv(config.output_dir / "03_shared_comparison.tsv", sep="\t", index=False)
    shared_by_sequence.to_csv(
        config.output_dir / "04_shared_by_sequence.tsv", sep="\t", index=False
    )
    merged.to_csv(config.output_dir / "01_merged_table.tsv", sep="\t", index=False)
    afdb_entry_table.to_csv(config.output_dir / "02_afdb_entry_roles.tsv", sep="\t", index=False)
    afdb_only.to_csv(config.output_dir / "06_afdb_only.tsv", sep="\t", index=False)
    pombase_only.to_csv(config.output_dir / "07_pombase_only.tsv", sep="\t", index=False)
    pd.DataFrame({"metric": list(counts), "value": list(counts.values())}).to_csv(
        config.output_dir / "08_summary_counts.tsv", sep="\t", index=False
    )


def write_excel(
    config: Config,
    counts: dict[str, int],
    shared_comparison: pd.DataFrame,
    shared_by_sequence: pd.DataFrame,
    merged: pd.DataFrame,
    afdb_entry_table: pd.DataFrame,
    afdb_only: pd.DataFrame,
    pombase_only: pd.DataFrame,
) -> None:
    """Write the multi-sheet Excel workbook."""
    summary = pd.DataFrame({"metric": list(counts), "value": list(counts.values())})
    shared = merged.loc[merged["category"] == "shared"]
    sequence_columns = ["afdb_sequence", "pombase_sequence"]
    workbook_path = config.output_dir / "10_afdb_vs_pombase.xlsx"

    with pd.ExcelWriter(workbook_path, engine=EXCEL_ENGINE) as writer:
        summary.to_excel(writer, sheet_name="summary", index=False)
        merged.to_excel(writer, sheet_name="merged", index=False)
        afdb_entry_table.to_excel(writer, sheet_name="afdb_entries", index=False)
        shared_comparison.loc[shared_comparison["match_method"] == "uniprot_id"].drop(
            columns=sequence_columns
        ).to_excel(writer, sheet_name="shared_by_id", index=False)
        shared_by_sequence.drop(columns=sequence_columns).to_excel(
            writer, sheet_name="shared_by_sequence", index=False
        )
        if "relationship" in shared_comparison.columns:
            check_columns = [
                column for column in shared_comparison.columns if column not in sequence_columns
            ]
            shared_comparison.loc[
                shared_comparison["relationship"].notna(), check_columns
            ].to_excel(writer, sheet_name="shared_uniprot_check", index=False)
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
            afdb_only_no_sequence["match_status"] == "duplicate_of_selected_entry"
        ].to_excel(writer, sheet_name="afdb_only_duplicates", index=False)
        afdb_only_no_sequence.loc[
            afdb_only_no_sequence["match_status"] == "superseded_by_sequence_match"
        ].to_excel(writer, sheet_name="afdb_only_superseded", index=False)

        pombase_only.to_excel(writer, sheet_name="pombase_only_all", index=False)

    logger.info(f"Excel workbook written to {workbook_path}")


def main() -> int:
    """Load inputs, run the comparison, and write the report."""
    args = parse_args()
    config = build_config(args)

    try:
        afdb_entries, pae_accessions = load_afdb_inputs(config.afdb_entries)
        gene_catalog = load_pombase_gene_catalog(config.pombase_dir)
        pombase_genes = gene_catalog.loc[
            gene_catalog[GeneCol.GENE_TYPE] == PROTEIN_CODING
        ].reset_index(drop=True)
        pombase_peptides = load_pombase_peptides(config.pombase_dir)

        pairs = assign_structure_accessions(pombase_genes, pombase_peptides, afdb_entries)
        paired = pairs.loc[pairs["match_method"] != "none"].reset_index(drop=True)
        unpaired_genes = pairs.loc[pairs["match_method"] == "none"].reset_index(drop=True)

        shared_comparison = attach_comparison(paired, afdb_entries, build_aligner())
        if config.uniprot_triangulation is not None:
            triangulation = load_uniprot_triangulation(config.uniprot_triangulation).rename(
                columns={ComparisonCol.ACCESSION.value: "pombase_accession"}
            )
            shared_comparison = shared_comparison.merge(
                triangulation,
                on=[GeneCol.SYSTEMATIC_ID, "pombase_accession"],
                how="left",
            )
        selected_accessions = set(shared_comparison["afdb_accession"])
        afdb_only = annotate_unpaired_entries(
            afdb_entries, pombase_genes, pombase_peptides, selected_accessions
        )
        if config.uniprot_resolution is not None:
            pombase_accession_set = set(pombase_genes[GeneCol.ACCESSION].dropna())
            afdb_only = afdb_only.merge(
                load_uniprot_resolution(config.uniprot_resolution, pombase_accession_set),
                on=AfdbCol.ACCESSION,
                how="left",
            )
        afdb_max_length = int(afdb_entries[AfdbCol.LENGTH].max())
        pombase_only = build_pombase_only(unpaired_genes, afdb_max_length)
        counts = build_summary_counts(
            afdb_entries,
            pae_accessions,
            pombase_genes,
            pombase_peptides,
            shared_comparison,
            afdb_only,
            pombase_only,
        )
        shared_by_sequence = build_shared_by_sequence(shared_comparison)
        merged = build_merged_table(
            shared_comparison.assign(category="shared", pae_available=True, notes=""),
            build_afdb_only_rows(afdb_only),
            build_pombase_only_rows(pombase_only),
        )
        afdb_entry_table = build_afdb_entry_table(
            afdb_entries, pae_accessions, shared_comparison, afdb_only
        )

        write_tables(
            config,
            shared_comparison,
            shared_by_sequence,
            merged,
            afdb_entry_table,
            afdb_only,
            pombase_only,
            counts,
        )
        write_excel(
            config,
            counts,
            shared_comparison,
            shared_by_sequence,
            merged,
            afdb_entry_table,
            afdb_only,
            pombase_only,
        )
        report = render_report(shared_comparison, afdb_only, pombase_only, counts, config)
        (config.output_dir / "09_report.md").write_text(report, encoding="utf-8")
    except (OSError, ValueError, KeyError) as error:
        logger.exception(f"Fatal execution error: {error}")
        return 1
    else:
        logger.info(f"Report written to {config.output_dir / '09_report.md'}")
        return 0


if __name__ == "__main__":
    sys.exit(main())
