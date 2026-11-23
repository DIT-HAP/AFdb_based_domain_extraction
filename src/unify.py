"""
Unified segment table.
======================

Projects the three domain sources onto a single long-format schema:

- ``PAE``  — AlphaFold PAE clustering merged with pLDDT confidence.
- ``TED``  — TED domains for the PomBase proteome.
- ``ECOD`` — ECOD domains for the PomBase proteome.

TED and ECOD are filtered to the proteins selected by the registry, so the
table never carries external segments for proteins the registry did not model.
All coordinates are converted to 1-based inclusive. Multi-interval rows are
exploded into one row per interval, and ``segment_part`` records the position
within the original domain.

Input
-----
- Annotated PAE table, PomBase TED table, PomBase ECOD table, registry.

Output
------
- ``results/11_segments.tsv``.
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd

from src.files import read_tsv, require_columns, write_tsv

# =============================================================================
# 2. CONSTANTS & COLUMN CONTRACTS
# =============================================================================
UNIFIED_COLUMNS = [
    "dataset",
    "uniprot_id",
    "systematic_id",
    "gene_name",
    "protein_length",
    "source",
    "source_version",
    "segment_id",
    "segment_part",
    "start",
    "end",
    "length",
    "score",
    "score_type",
    "fraction_gt_70",
    "segment_class",
    "family",
    "location",
]
PAE_REQUIRED = [
    "Systematic ID",
    "uniprot_id",
    "protein_length",
    "domain_id",
    "domain_residues",
    "mean_pLDDT",
    "fraction_gt_70",
    "domain_type",
    "location",
]
TED_REQUIRED = ["ted_id", "uniprot_id", "consensus_level", "chopping", "plddt", "cath_label"]
ECOD_REQUIRED = ["Protein", "Column2", "Range", "ID", "uniprot_id"]
REGISTRY_REQUIRED = [
    "gene_systematic_id",
    "uniprot_id",
    "gene_name",
    "modeled_length",
    "ted_ecod_coordinate_ok",
    "in_analysis",
]
SOURCE_ORDER = ["PAE", "TED", "ECOD"]


# =============================================================================
# 3. FUNCTIONAL CORE
# =============================================================================
def expand_intervals(
    frame: pd.DataFrame,
    interval_column: str,
    separator: str,
    group_columns: tuple[str, ...] = ("uniprot_id", "segment_id"),
) -> pd.DataFrame:
    """Explode a multi-interval column into one 1-based row per interval."""
    expanded = frame.assign(_interval=frame[interval_column].str.split(separator)).explode(
        "_interval"
    )
    bounds = expanded["_interval"].str.split("-", expand=True)
    expanded["start"] = bounds[0].astype(int)
    expanded["end"] = bounds[1].astype(int)
    expanded["length"] = expanded["end"] - expanded["start"] + 1
    expanded["segment_part"] = expanded.groupby(list(group_columns), sort=False).cumcount() + 1
    return expanded.drop(columns=["_interval", interval_column])


def load_pae_segments(path: Path, dataset: str, source_version: str) -> pd.DataFrame:
    """Load the annotated PAE table into the unified schema."""
    frame = read_tsv(path)
    require_columns(frame, PAE_REQUIRED, path)
    frame = frame.rename(
        columns={
            "Systematic ID": "systematic_id",
            "domain_id": "segment_id",
            "mean_pLDDT": "score",
            "domain_type": "segment_class",
        }
    )
    frame["segment_id"] = frame["segment_id"].astype(str)
    frame["dataset"] = dataset
    frame["source"] = "PAE"
    frame["source_version"] = source_version
    frame["score_type"] = "mean_pLDDT"
    frame["family"] = ""
    frame["gene_name"] = pd.NA
    return expand_intervals(
        frame,
        "domain_residues",
        ",",
        group_columns=("systematic_id", "uniprot_id", "segment_id"),
    )[UNIFIED_COLUMNS]


def load_ted_segments(path: Path, dataset: str, source_version: str) -> pd.DataFrame:
    """Load the PomBase TED table into the unified schema."""
    frame = read_tsv(path)
    require_columns(frame, TED_REQUIRED, path)
    frame = frame.rename(
        columns={
            "ted_id": "segment_id",
            "plddt": "score",
            "consensus_level": "segment_class",
            "cath_label": "family",
        }
    )
    frame["dataset"] = dataset
    frame["source"] = "TED"
    frame["source_version"] = source_version
    frame["score_type"] = "TED_plddt"
    frame["fraction_gt_70"] = np.nan
    frame["location"] = ""
    frame["systematic_id"] = pd.NA
    frame["gene_name"] = pd.NA
    frame["protein_length"] = pd.NA
    return expand_intervals(frame, "chopping", "_")[UNIFIED_COLUMNS]


def load_ecod_segments(path: Path, dataset: str, source_version: str) -> pd.DataFrame:
    """Load the PomBase ECOD table into the unified schema."""
    frame = read_tsv(path)
    require_columns(frame, ECOD_REQUIRED, path)
    frame["segment_id"] = frame["Protein"] + "_" + frame["Column2"]
    frame = frame.rename(columns={"Column2": "segment_class", "ID": "family"})
    frame["dataset"] = dataset
    frame["source"] = "ECOD"
    frame["source_version"] = source_version
    frame["score"] = np.nan
    frame["score_type"] = "NA"
    frame["fraction_gt_70"] = np.nan
    frame["location"] = ""
    frame["systematic_id"] = pd.NA
    frame["gene_name"] = pd.NA
    frame["protein_length"] = pd.NA
    return expand_intervals(frame, "Range", ",")[UNIFIED_COLUMNS]


def attach_pombase_annotation(
    frame: pd.DataFrame,
    genes: pd.DataFrame,
    stats: pd.DataFrame,
) -> pd.DataFrame:
    """Fill ids, gene names, and protein lengths from PomBase, falling back to the systematic id."""
    gene_map = genes[["uniprot_id", "gene_systematic_id", "gene_name"]].drop_duplicates(
        "uniprot_id"
    )
    length_map = stats[["Systematic_ID", "Residues"]].rename(
        columns={"Systematic_ID": "gene_systematic_id", "Residues": "pombase_length"}
    )
    gene_map = gene_map.merge(length_map, on="gene_systematic_id", how="left")

    result = frame.merge(gene_map, on="uniprot_id", how="left", suffixes=("", "_pombase"))
    result["systematic_id"] = result["systematic_id"].combine_first(result["gene_systematic_id"])
    names_by_gene = genes.drop_duplicates("gene_systematic_id").set_index("gene_systematic_id")
    result["gene_name"] = (
        result["gene_name"]
        .combine_first(result["systematic_id"].map(names_by_gene["gene_name"]))
        .combine_first(result["gene_name_pombase"])
    )
    blank_name = result["gene_name"].fillna("").eq("")
    result["gene_name"] = result["gene_name"].mask(blank_name, result["systematic_id"])
    result["protein_length"] = result["protein_length"].combine_first(result["pombase_length"])
    result = result.drop(columns=["gene_systematic_id", "gene_name_pombase", "pombase_length"])
    result["protein_length"] = result["protein_length"].astype("Int64")
    return result


def _sort_segments(frame: pd.DataFrame) -> pd.DataFrame:
    """Sort segments by dataset, protein, source, and coordinate."""
    result = frame.copy()
    result["source"] = pd.Categorical(result["source"], categories=SOURCE_ORDER, ordered=True)
    return result.sort_values(
        ["dataset", "uniprot_id", "source", "start", "segment_part"],
        kind="stable",
        ignore_index=True,
    )


# =============================================================================
# 4. IMPERATIVE SHELL
# =============================================================================
def unify_segments_from_registry(
    pae_path: Path,
    ted_path: Path,
    ecod_path: Path,
    registry_path: Path,
    output_path: Path,
    *,
    dataset: str,
    pae_source_version: str,
    ted_version: str,
    ecod_version: str,
) -> None:
    """Build the unified segment table from registry-selected domains."""
    pae = load_pae_segments(pae_path, dataset, pae_source_version)
    registry = read_tsv(registry_path)
    require_columns(registry, REGISTRY_REQUIRED, registry_path)

    coordinate_ok = registry.loc[registry["in_analysis"] & registry["ted_ecod_coordinate_ok"]]
    dataset_proteins = set(coordinate_ok["uniprot_id"].astype(str))
    ted = load_ted_segments(ted_path, dataset, ted_version)
    ecod = load_ecod_segments(ecod_path, dataset, ecod_version)
    ted = ted.loc[ted["uniprot_id"].astype(str).isin(dataset_proteins)]
    ecod = ecod.loc[ecod["uniprot_id"].astype(str).isin(dataset_proteins)]
    combined = pd.concat([pae, ted, ecod], ignore_index=True)

    genes = registry[["uniprot_id", "gene_systematic_id", "gene_name"]]
    stats = registry[["gene_systematic_id", "modeled_length"]].rename(
        columns={"gene_systematic_id": "Systematic_ID", "modeled_length": "Residues"}
    )
    annotated = attach_pombase_annotation(combined, genes, stats)
    write_tsv(_sort_segments(annotated)[UNIFIED_COLUMNS], output_path)
