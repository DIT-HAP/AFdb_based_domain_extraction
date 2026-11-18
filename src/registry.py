"""
Gene structure registry.
========================

Selects exactly one structure and one PAE matrix per PomBase protein-coding
gene for the pinned PomBase 2026-06-01 / AFdb v6 comparison.

Selection rules, in priority order:

1. The AFdb v6 entry already selected by the sequence-first merge (merged table
   ``afdb_accession``) is used when its modelled sequence equals the PomBase
   peptide. This covers both the gene's own UniProt entry and sequence-identical
   alternate accessions.
2. Remaining genes use the local vDulab prediction (2024-03-05) when its input
   sequence matches the PomBase peptide.
3. The three long proteins without vDulab use the local AF3 prediction
   (2024-11-01).
4. Everything else is excluded with an explicit reason; AFdb entries without a
   selected gene stay in the table but never enter the downstream analysis.

Input
-----
- Comparison merged table with the selected entries, both sequences, and the
  UniProt relationship column.
- AFdb v4 CIF files, vDulab structures/PAE, AF3 structures/PAE.

Output
------
- ``gene_structure_registry.tsv`` with one row per PomBase coding gene plus
  the AFdb-only entries, ``registry_summary.tsv`` with selection counts, and
  ``registry_stats.tsv`` with downstream counts.
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from typing import TextIO

from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB.MMCIFParser import MMCIFParser
from loguru import logger
import numpy as np
import pandas as pd

from src.afdb import read_cif_header_sequence
from src.files import AppConfig, read_tsv, require_columns, write_tsv

# =============================================================================
# 2. CONSTANTS & COLUMN CONTRACTS
# =============================================================================
CIF_NAME_PATTERN = re.compile(r"^AF-(?P<accession>.+)-F\d+-model_v[^.]+\.cif(?:\.gz)?$")
AFDB_PRIMARY = "AFDB_v6_primary"
AFDB_SEQUENCE_MATCH = "AFDB_v6_sequence_match"
AFDB_ALTERNATE = "AFDB_v6_alternate"
VDULAB_SOURCE = "vDulab"
AF3_SOURCE = "AF3"
NO_SOURCE = "none"
NONIDENTICAL_CATEGORIES = ("substitution_only", "indel_only", "both")
EXCLUSION_REASONS = {
    "afdb_only_not_analysed",
    "excluded_partial_structure",
    "excluded_no_structure",
    "excluded_prediction_sequence_mismatch",
    "excluded_no_valid_prediction",
}
POMBASE_ONLY_NO_STRUCTURE = {
    "SPBC460.01c",
    "SPBC460.02c",
    "SPBC460.03",
    "SPBC460.04c",
    "SPBC460.05",
}
PARTIAL_STRUCTURE_GENE = "SPAC23G3.02c"
THIRD_VERSION_GENE = "SPCC569.03"

MERGED_REQUIRED = [
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
    "afdb_sequence",
    "afdb_length",
    "pombase_length",
    "pombase_sequence",
    "length_delta",
    "sequence_category",
    "relationship",
    "n_substitutions",
    "n_afdb_extra",
    "n_pombase_extra",
    "identity_pct",
    "duplicate_afdb_accessions",
    "notes",
]

REGISTRY_COLUMNS = [
    "gene_systematic_id",
    "gene_name",
    "gene_product",
    "gene_type",
    "uniprot_id",
    "pombase_accession",
    "afdb_accession",
    "afdb_description",
    "afdb_length",
    "afdb_sequence",
    "pombase_length",
    "pombase_sequence",
    "category",
    "match_method",
    "sequence_category",
    "relationship",
    "length_delta",
    "n_substitutions",
    "n_afdb_extra",
    "n_pombase_extra",
    "identity_pct",
    "duplicate_afdb_accessions",
    "selection_reason",
    "in_analysis",
    "exclusion_reason",
    "structure_source",
    "structure_id",
    "structure_path",
    "pae_path",
    "modeled_length",
    "structure_sequence",
    "coordinate_basis",
    "plddt_granularity",
    "ted_ecod_coordinate_ok",
    "notes",
]


# =============================================================================
# 3. CONFIGURATION
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class RegistryPaths:
    """Immutable paths and release labels for a registry build."""

    merged_table: Path
    afdb_v4_dir: Path
    vdulab_structure_dir: Path
    vdulab_fasta_dir: Path
    af3_structure_dir: Path
    output: Path
    summary_output: Path
    stats_output: Path


# =============================================================================
# 4. BOUNDARY LOADERS
# =============================================================================
def open_text(path: Path) -> TextIO:
    """Open a plain or gzip-compressed text file."""
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def _header_sequence(path: Path) -> str:
    """Read the modelled sequence from an AFdb CIF header when present."""
    with open_text(path) as handle:
        return read_cif_header_sequence(handle)[1]


def _atom_sequence(path: Path) -> str:
    """Read a CIF sequence from its CA atoms when the header omits it."""
    with open_text(path) as handle:
        structure = MMCIFParser(QUIET=True).get_structure(path.stem, handle)
    return "".join(
        protein_letters_3to1.get(residue.get_resname(), "X")
        for residue in structure.get_residues()
        if residue.has_id("CA")
    )


def read_cif_sequence(path: Path) -> str:
    """Read the modelled sequence from a CIF header or its CA atoms."""
    return _header_sequence(path) or _atom_sequence(path)


def load_cif_sequences(directory: Path) -> dict[str, str]:
    """Map every model accession in a directory to its modelled sequence."""
    sequences: dict[str, str] = {}
    for path in sorted(directory.glob("AF-*-model_v*.cif*")):
        match = CIF_NAME_PATTERN.match(path.name)
        if match is None:
            raise ValueError(f"Unrecognised CIF filename: {path.name}")
        sequences[match.group("accession")] = read_cif_sequence(path)
    if not sequences:
        raise FileNotFoundError(f"No CIF structures found in {directory}")
    return sequences


def load_gene_fastas(directory: Path) -> dict[str, str]:
    """Map gene ids to the single sequence of their prediction FASTA."""
    sequences: dict[str, str] = {}
    for path in sorted(directory.rglob("*.fasta")):
        records = path.read_text().split(">")[1:]
        if len(records) != 1:
            raise ValueError(f"Expected one record in {path.name}, found {len(records)}")
        sequences[path.stem] = records[0].split("\n", 1)[1].replace("\n", "").rstrip("*")
    return sequences


def load_structure_ids(directory: Path, suffix: str, model_tag: str) -> list[str]:
    """List the model identifiers of one custom prediction directory."""
    pattern = f"AF-*-F1-model_{model_tag}.{suffix}"
    return [
        path.name.removeprefix("AF-").split("-F1-")[0] for path in sorted(directory.glob(pattern))
    ]


# =============================================================================
# 5. FUNCTIONAL CORE
# =============================================================================
def _match_prediction(
    merged: pd.DataFrame,
    stems: list[str],
    prediction_sequences: dict[str, str],
    *,
    sequence_key: str = "gene_systematic_id",
) -> pd.DataFrame:
    """Return genes whose prediction sequence equals the PomBase peptide.

    ``sequence_key`` is the gene id for vDulab FASTAs and the UniProt accession
    for AF3 outputs.
    """
    available = merged[["gene_systematic_id", "uniprot_id", "pombase_sequence"]].copy()
    available["structure_id"] = available["uniprot_id"].where(available["uniprot_id"].isin(stems))
    available["structure_sequence"] = available[sequence_key].map(prediction_sequences)
    return available.loc[
        available["structure_id"].notna()
        & (available["pombase_sequence"] == available["structure_sequence"]),
        ["gene_systematic_id", "structure_id", "structure_sequence"],
    ]


def _selection_reason(
    merged: pd.DataFrame,
    vdulab: pd.DataFrame,
    af3: pd.DataFrame,
) -> np.ndarray:
    """Assign one selection reason per row with named boolean masks."""
    is_afdb_only = merged["category"] == "afdb_only"
    is_identical = merged["sequence_category"] == "identical"
    is_primary = is_identical & merged["afdb_accession"].eq(merged["uniprot_id"])
    is_alternate = is_identical & ~is_primary
    has_vdulab = merged["gene_systematic_id"].isin(vdulab["gene_systematic_id"])
    has_af3 = merged["gene_systematic_id"].isin(af3["gene_systematic_id"])
    is_partial = merged["gene_systematic_id"] == PARTIAL_STRUCTURE_GENE
    has_no_structure = merged["gene_systematic_id"].isin(POMBASE_ONLY_NO_STRUCTURE)
    is_third_version = merged["gene_systematic_id"] == THIRD_VERSION_GENE

    conditions = [
        is_afdb_only,
        is_primary,
        is_alternate,
        has_vdulab,
        has_af3,
        is_partial,
        has_no_structure,
        is_third_version,
    ]
    choices = [
        "afdb_only_not_analysed",
        AFDB_PRIMARY,
        AFDB_ALTERNATE,
        VDULAB_SOURCE,
        AF3_SOURCE,
        "excluded_partial_structure",
        "excluded_no_structure",
        "excluded_prediction_sequence_mismatch",
    ]
    return np.select(conditions, choices, default="excluded_no_valid_prediction")


def _attach_structure(
    registry: pd.DataFrame,
    vdulab: pd.DataFrame,
    af3: pd.DataFrame,
    v4_sequences: dict[str, str],
) -> pd.DataFrame:
    """Attach the selected structure id and modelled sequence to every row."""
    is_identical = registry["sequence_category"].eq("identical")
    afdb_sources = pd.DataFrame(
        {
            "gene_systematic_id": registry["gene_systematic_id"],
            "structure_id": registry["afdb_accession"],
            "structure_sequence": registry["afdb_sequence"],
        }
    )
    afdb_sources = afdb_sources.loc[is_identical & afdb_sources["structure_id"].notna()]
    sources = pd.concat([afdb_sources, vdulab, af3], ignore_index=True)
    sources = sources.drop_duplicates("gene_systematic_id", keep="first")

    result = registry.merge(sources, on="gene_systematic_id", how="left")
    is_primary = result["selection_reason"].eq(AFDB_PRIMARY)
    has_pombase_accession = result["pombase_accession"].notna()
    is_sequence_match = is_primary & (
        ~has_pombase_accession | result["pombase_accession"].ne(result["uniprot_id"])
    )
    result["structure_source"] = np.select(
        [
            is_primary & ~is_sequence_match,
            is_sequence_match,
            result["selection_reason"].eq(AFDB_ALTERNATE),
            result["selection_reason"].eq(VDULAB_SOURCE),
            result["selection_reason"].eq(AF3_SOURCE),
        ],
        [AFDB_PRIMARY, AFDB_SEQUENCE_MATCH, AFDB_ALTERNATE, VDULAB_SOURCE, AF3_SOURCE],
        default=NO_SOURCE,
    )
    result["modeled_length"] = result["structure_sequence"].str.len().astype("Int64")
    result["coordinate_basis"] = np.select(
        [
            result["structure_source"].isin([AFDB_PRIMARY, AFDB_SEQUENCE_MATCH]),
            result["structure_source"].eq(AFDB_ALTERNATE),
            result["structure_source"].isin([VDULAB_SOURCE, AF3_SOURCE]),
        ],
        ["UniProt", "PomBase", "PomBase"],
        default="",
    )
    result["plddt_granularity"] = np.where(
        result["structure_source"].eq(AF3_SOURCE), "atom_ca", "residue"
    )
    no_source = result["structure_source"].eq(NO_SOURCE)
    result.loc[no_source, ["coordinate_basis", "plddt_granularity"]] = None
    result["in_analysis"] = ~result["selection_reason"].isin(EXCLUSION_REASONS)
    result["ted_ecod_coordinate_ok"] = (
        result["structure_sequence"].eq(result["uniprot_id"].map(v4_sequences))
        & result["in_analysis"]
    )
    result["exclusion_reason"] = np.where(result["in_analysis"], "", result["selection_reason"])
    return result


def _attach_paths(registry: pd.DataFrame) -> pd.DataFrame:
    """Build config-resolvable structure and PAE paths for each selected source."""
    structure_id = registry["structure_id"]
    source = registry["structure_source"]
    is_afdb = source.isin([AFDB_PRIMARY, AFDB_SEQUENCE_MATCH, AFDB_ALTERNATE])
    is_vdulab = source.eq(VDULAB_SOURCE)
    is_af3 = source.eq(AF3_SOURCE)
    registry["structure_path"] = np.select(
        [is_afdb, is_vdulab, is_af3],
        [
            "reference:UP000002485_284812_SCHPO_v6/AF-" + structure_id + "-F1-model_v6.cif.gz",
            "raw:20240305_predicted_results/vDulab/AF-" + structure_id + "-F1-model_vDulab.pdb.gz",
            "raw:20241101_AF3_output/PDB/AF-" + structure_id + "-F1-model_vAF3.cif",
        ],
        default="",
    )
    registry["pae_path"] = np.select(
        [is_afdb, is_vdulab, is_af3],
        [
            "reference:PAE/UP000002485_284812_SCHPO_v6/AF-"
            + structure_id
            + "-F1-predicted_aligned_error_v6.json",
            "raw:20240305_predicted_results/PAE_files/AF-"
            + structure_id
            + "-F1-predicted_aligned_error_vDulab.json",
            "raw:20241101_AF3_output/PAE_files/AF-"
            + structure_id
            + "-F1-predicted_aligned_error_vAF3.json",
        ],
        default="",
    )
    return registry


def build_registry_stats(registry: pd.DataFrame) -> pd.DataFrame:
    """Summarise the registry-derived downstream counts."""
    in_analysis = registry["in_analysis"]
    source = registry["structure_source"]
    is_nonidentical = registry["sequence_category"].isin(NONIDENTICAL_CATEGORIES)
    metrics = {
        "in_analysis": int(in_analysis.sum()),
        "nonidentical_shared": int((is_nonidentical & registry["category"].eq("shared")).sum()),
        "vdulab_total": int((source.eq(VDULAB_SOURCE) & in_analysis).sum()),
        "vdulab_nonidentical": int(
            (source.eq(VDULAB_SOURCE) & in_analysis & is_nonidentical).sum()
        ),
        "vdulab_long": int(
            (
                source.eq(VDULAB_SOURCE)
                & in_analysis
                & registry["sequence_category"].eq("not_compared")
            ).sum()
        ),
        "af3_total": int((source.eq(AF3_SOURCE) & in_analysis).sum()),
        "excluded_prediction_sequence_mismatch": int(
            registry["selection_reason"].eq("excluded_prediction_sequence_mismatch").sum()
        ),
        "ted_ecod_coordinate_ok": int((in_analysis & registry["ted_ecod_coordinate_ok"]).sum()),
    }
    return pd.DataFrame({"metric": list(metrics), "value": list(metrics.values())})


# =============================================================================
# 6. IMPERATIVE SHELL
# =============================================================================
def _missing_paths(specs: pd.Series, app_config: AppConfig) -> list[str]:
    """Return the path specifications that do not resolve to an existing file."""
    return [spec for spec in specs if not app_config.resolve(spec).exists()]


def _validate_selection(registry: pd.DataFrame, app_config: AppConfig) -> None:
    """Fail fast when a selected row lacks a readable structure or PAE file."""
    selected = registry.loc[registry["in_analysis"]]
    missing_structure = _missing_paths(selected["structure_path"], app_config)
    missing_pae = _missing_paths(selected["pae_path"], app_config)
    if missing_structure or missing_pae:
        raise FileNotFoundError(
            f"Missing selected files: structures={missing_structure[:5]}, PAE={missing_pae[:5]}"
        )


def build_registry(config: RegistryPaths, app_config: AppConfig) -> None:
    """Build and write the gene structure registry, summary, and stats."""
    registry = read_tsv(config.merged_table)
    require_columns(registry, MERGED_REQUIRED, config.merged_table)

    v4_sequences = load_cif_sequences(config.afdb_v4_dir)
    af3_sequences = load_cif_sequences(config.af3_structure_dir)
    vdulab_sequences = load_gene_fastas(config.vdulab_fasta_dir)
    vdulab_stems = load_structure_ids(config.vdulab_structure_dir, "pdb.gz", VDULAB_SOURCE)
    af3_stems = load_structure_ids(config.af3_structure_dir, "cif", "vAF3")

    vdulab = _match_prediction(registry, vdulab_stems, vdulab_sequences)
    af3 = _match_prediction(registry, af3_stems, af3_sequences, sequence_key="uniprot_id")

    registry["selection_reason"] = _selection_reason(registry, vdulab, af3)
    registry = _attach_structure(registry, vdulab, af3, v4_sequences)
    registry = _attach_paths(registry)
    _validate_selection(registry, app_config)

    output = registry[REGISTRY_COLUMNS].sort_values("gene_systematic_id", ignore_index=True)
    write_tsv(output, config.output)

    summary = (
        output.groupby(["selection_reason", "structure_source", "in_analysis"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False, ignore_index=True)
    )
    write_tsv(summary, config.summary_output)
    write_tsv(build_registry_stats(output), config.stats_output)
    logger.info(
        f"Registry written to {config.output}: {int(output['in_analysis'].sum()):,} "
        f"of {len(output):,} rows in analysis"
    )
