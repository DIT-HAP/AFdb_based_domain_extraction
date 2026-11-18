"""
PAE domain merging and annotation.
==================================

Implements the three-step domain refinement that defined the legacy results:

1. Merge sub-domains whose mean pLDDT is below the threshold into disordered
   regions.
2. Merge confident sub-domains of one cluster when they are separated by fewer
   than ``merge_gap`` residues.
3. Remove merged regions that overlap with a confident contiguous domain.

The same module also implements the independent 3D-AF-Surfer style confident
domain caller.

Input
-----
- PAE domain table, pLDDT table, and a PomBase gene table.

Output
------
- Reformatted domain table and its annotated variant.
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd

from src.files import PipelineParameters, read_tsv, require_columns, write_tsv
from src.plddt import check_overlap, domain_stats, merge_intervals

# =============================================================================
# 2. CONSTANTS & COLUMN CONTRACTS
# =============================================================================
PAE_REQUIRED = ["uniprot_id", "version", "domain_id", "domain_residues", "protein_length"]
PLDDT_REQUIRED = ["uniprot_id", "protein_length", "pLDDT"]
GENE_REQUIRED = ["gene_systematic_id", "uniprot_id"]
REGISTRY_REQUIRED = ["gene_systematic_id", "uniprot_id"]
REGISTRY_PAE_REQUIRED = [
    "gene_systematic_id",
    "version",
    "domain_id",
    "domain_residues",
    "protein_length",
]
REGISTRY_PLDDT_REQUIRED = ["gene_systematic_id", "protein_length", "pLDDT"]
STAT_COLUMNS = ["domain_length", "mean_pLDDT", "fraction_gt_70"]
REFORMATED_COLUMNS = [
    "Systematic ID",
    "uniprot_id",
    "version",
    "domain_id",
    "domain_residues",
    "protein_length",
    *STAT_COLUMNS,
]
CONFIDENT_COLUMNS = [
    "Systematic ID",
    "uniprot_id",
    "domain_id",
    "domain_residues",
    "domain_type",
    "protein_length",
    *STAT_COLUMNS,
]


# =============================================================================
# 3. FUNCTIONAL CORE (Pure Transformations)
# =============================================================================
def merge_domain_group(group: pd.DataFrame, parameters: PipelineParameters) -> list[str]:
    """Apply the legacy pLDDT merging rules to one protein."""
    if len(group) == 1:
        return group["domain_residues"].tolist()

    is_contiguous = group["domain_residues"].str.contains(",") & (
        group["mean_pLDDT"] >= parameters.plddt_threshold
    )
    contiguous_count = int(is_contiguous.sum())
    if contiguous_count == 0:
        return _merge_non_contiguous(group, parameters)
    if contiguous_count == 1:
        return _merge_single_contiguous(group, is_contiguous, parameters)
    return _merge_multi_contiguous(group, is_contiguous, parameters)


def _disordered_intervals(group: pd.DataFrame, threshold: float) -> list[str]:
    """Flatten all intervals of sub-domains below the confidence threshold."""
    is_disordered = group["mean_pLDDT"] < threshold
    return group.loc[is_disordered, "domain_residues"].str.split(",").explode().tolist()


def _merge_non_contiguous(group: pd.DataFrame, parameters: PipelineParameters) -> list[str]:
    """Merge disordered sub-domains while keeping confident single intervals."""
    is_disordered = group["mean_pLDDT"] < parameters.plddt_threshold
    merged_disordered = merge_intervals(_disordered_intervals(group, parameters.plddt_threshold))
    confident = group.loc[~is_disordered, "domain_residues"].tolist()
    return merged_disordered + confident


def _merge_single_contiguous(
    group: pd.DataFrame,
    is_contiguous: pd.Series,
    parameters: PipelineParameters,
) -> list[str]:
    """Merge one confident contiguous cluster and drop regions it overlaps."""
    is_disordered = group["mean_pLDDT"] < parameters.plddt_threshold
    merged_disordered = merge_intervals(_disordered_intervals(group, parameters.plddt_threshold))
    confident_single = group.loc[~is_disordered & ~is_contiguous, "domain_residues"].tolist()

    contiguous_intervals = group.loc[is_contiguous, "domain_residues"].iloc[0].split(",")
    merged_contiguous = merge_intervals(contiguous_intervals, gap=parameters.merge_gap)

    filtered_disordered = [
        item for item in merged_disordered if not check_overlap(item, merged_contiguous)
    ]
    filtered_confident = [
        item for item in confident_single if not check_overlap(item, merged_contiguous)
    ]
    if len(merged_contiguous) == 1:
        return [*filtered_disordered, *filtered_confident, *merged_contiguous]
    return [*filtered_disordered, *filtered_confident, ",".join(merged_contiguous)]


def _merge_multi_contiguous(
    group: pd.DataFrame,
    is_contiguous: pd.Series,
    parameters: PipelineParameters,
) -> list[str]:
    """Merge several confident clusters, falling back when they overlap."""
    is_disordered = group["mean_pLDDT"] < parameters.plddt_threshold
    merged_disordered = merge_intervals(_disordered_intervals(group, parameters.plddt_threshold))
    confident_single = group.loc[~is_disordered & ~is_contiguous, "domain_residues"].tolist()
    contiguous_rows = group.loc[is_contiguous, "domain_residues"].tolist()

    merged_clusters = [
        merge_intervals(row.split(","), gap=parameters.merge_gap) for row in contiguous_rows
    ]
    flat_clusters = [interval for cluster in merged_clusters for interval in cluster]
    if _clusters_overlap(merged_clusters):
        return [*merged_disordered, *confident_single, *contiguous_rows]

    filtered_disordered = [
        item for item in merged_disordered if not check_overlap(item, flat_clusters)
    ]
    filtered_confident = [
        item for item in confident_single if not check_overlap(item, flat_clusters)
    ]
    return [
        *filtered_disordered,
        *filtered_confident,
        *(",".join(cluster) for cluster in merged_clusters),
    ]


def _clusters_overlap(merged_clusters: list[list[str]]) -> bool:
    """Return whether any merged cluster overlaps another cluster."""
    remaining = [interval for cluster in merged_clusters for interval in cluster]
    overlapping: list[str] = []
    for cluster in merged_clusters:
        remaining = [interval for interval in remaining if interval not in cluster]
        overlapping.extend(interval for interval in cluster if check_overlap(interval, remaining))
    return bool(set(overlapping))


def build_plddt_lookup(plddt: pd.DataFrame, id_column: str = "uniprot_id") -> dict[str, np.ndarray]:
    """Map each protein to its per-residue pLDDT array."""
    unique = plddt.drop_duplicates(id_column)
    return {
        protein_id: np.asarray(values.split(","), dtype=float)
        for protein_id, values in zip(unique[id_column], unique["pLDDT"], strict=True)
    }


def attach_stats(
    frame: pd.DataFrame,
    plddt_lookup: dict[str, np.ndarray],
    threshold: float,
    id_column: str = "uniprot_id",
) -> pd.DataFrame:
    """Attach length and pLDDT summary columns to a domain table."""
    stats = [
        domain_stats(interval, plddt_lookup[protein_id], threshold)
        for interval, protein_id in zip(frame["domain_residues"], frame[id_column], strict=True)
    ]
    result = frame.copy()
    result["domain_length"] = [row[0] for row in stats]
    result["mean_pLDDT"] = [row[1] for row in stats]
    result["fraction_gt_70"] = [row[2] for row in stats]
    return result


def _sort_and_number(frame: pd.DataFrame, id_column: str) -> pd.DataFrame:
    """Sort domains by genomic start and assign 1-based domain ids."""
    group_keys = [id_column, "version", "protein_length"]
    result = frame.copy()
    result["_start"] = result["domain_residues"].str.split("-").str[0].astype(int)
    result = result.sort_values([*group_keys, "_start"], kind="stable").drop(columns="_start")
    result["domain_id"] = result.groupby(group_keys, sort=False).cumcount() + 1
    return result


def _join_gene_annotation(
    frame: pd.DataFrame,
    genes: pd.DataFrame,
    *,
    fill_missing: bool,
) -> pd.DataFrame:
    """Attach PomBase systematic ids and gene names, optionally falling back to UniProt."""
    gene_map = genes[["gene_systematic_id", "uniprot_id"]].drop_duplicates("uniprot_id")
    result = gene_map.merge(frame, on="uniprot_id", how="right")
    if fill_missing:
        result["gene_systematic_id"] = np.where(
            result["gene_systematic_id"].isna(),
            result["uniprot_id"],
            result["gene_systematic_id"],
        )
    return result.rename(columns={"gene_systematic_id": "Systematic ID"})


def _join_registry_annotation(frame: pd.DataFrame, registry: pd.DataFrame) -> pd.DataFrame:
    """Attach the registry's gene id and primary accession."""
    annotation = registry[["gene_systematic_id", "uniprot_id"]]
    result = annotation.merge(frame, on="gene_systematic_id", how="right")
    return result.rename(columns={"gene_systematic_id": "Systematic ID"})


def build_reformated_frame(
    pae: pd.DataFrame,
    plddt: pd.DataFrame,
    parameters: PipelineParameters,
    *,
    id_column: str = "uniprot_id",
) -> pd.DataFrame:
    """Merge PAE clusters with pLDDT confidence into one reformatted frame."""
    lookup = build_plddt_lookup(plddt, id_column)
    merged = pae.merge(plddt[[id_column, "protein_length"]], on=[id_column, "protein_length"])
    merged = attach_stats(merged, lookup, parameters.plddt_threshold, id_column)
    exploded = _merge_groups(merged, parameters, id_column)
    numbered = _sort_and_number(exploded, id_column)
    return attach_stats(numbered, lookup, parameters.plddt_threshold, id_column)


# =============================================================================
# 4. IMPERATIVE SHELL (I/O Orchestration)
# =============================================================================
def _merge_groups(
    frame: pd.DataFrame, parameters: PipelineParameters, id_column: str
) -> pd.DataFrame:
    """Apply the merge rules per protein and return one row per merged domain."""
    group_keys = [id_column, "version", "protein_length"]
    records = []
    for indices in frame.groupby(group_keys, sort=False).indices.values():
        group = frame.iloc[indices]
        first = group.iloc[0]
        records.extend(
            {
                id_column: first[id_column],
                "version": first["version"],
                "protein_length": first["protein_length"],
                "domain_residues": interval,
            }
            for interval in merge_domain_group(group, parameters)
        )
    return pd.DataFrame.from_records(records)


def build_reformated_table(
    pae_path: Path,
    plddt_path: Path,
    gene_path: Path,
    parameters: PipelineParameters,
    output_path: Path,
) -> None:
    """Merge PAE clusters with pLDDT confidence and write the reformatted table."""
    pae = read_tsv(pae_path)
    plddt = read_tsv(plddt_path)
    genes = read_tsv(gene_path)
    require_columns(pae, PAE_REQUIRED, pae_path)
    require_columns(plddt, PLDDT_REQUIRED, plddt_path)
    require_columns(genes, GENE_REQUIRED, gene_path)

    numbered = build_reformated_frame(pae, plddt, parameters)
    annotated = _join_gene_annotation(numbered, genes, fill_missing=True)
    write_tsv(annotated[REFORMATED_COLUMNS], output_path)


def build_reformated_table_from_registry(
    pae_path: Path,
    plddt_path: Path,
    registry_path: Path,
    parameters: PipelineParameters,
    output_path: Path,
) -> None:
    """Merge registry-selected PAE clusters with pLDDT and attach gene metadata."""
    pae = read_tsv(pae_path)
    plddt = read_tsv(plddt_path)
    registry = read_tsv(registry_path)
    require_columns(pae, REGISTRY_PAE_REQUIRED, pae_path)
    require_columns(plddt, REGISTRY_PLDDT_REQUIRED, plddt_path)
    require_columns(registry, REGISTRY_REQUIRED, registry_path)

    numbered = build_reformated_frame(pae, plddt, parameters, id_column="gene_systematic_id")
    annotated = _join_registry_annotation(numbered, registry)
    write_tsv(annotated[REFORMATED_COLUMNS], output_path)


def assign_type(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Assign N/C-terminal location and domain/disordered type."""
    first_interval = frame["domain_residues"].str.split(",").str[0]
    last_interval = frame["domain_residues"].str.split(",").str[-1]
    start = first_interval.str.split("-").str[0].astype(int)
    end = last_interval.str.split("-").str[1].astype(int)
    spans_whole = frame["domain_length"] == frame["protein_length"]
    at_n_terminus = (start == 1) & (end != frame["protein_length"])
    at_c_terminus = (start != 1) & (end == frame["protein_length"])

    location = np.select(
        [spans_whole, at_n_terminus, at_c_terminus],
        ["Whole", "N-terminal", "C-terminal"],
        default="middle",
    )
    domain_type = np.where(frame["mean_pLDDT"] > threshold, "domain", "disordered")
    return frame.assign(location=location, domain_type=domain_type)


def annotate_table(input_path: Path, parameters: PipelineParameters, output_path: Path) -> None:
    """Read the reformatted table, assign types, and write the annotated table."""
    frame = read_tsv(input_path)
    require_columns(frame, [*REFORMATED_COLUMNS], input_path)
    annotated = assign_type(frame, parameters.plddt_threshold)
    write_tsv(annotated, output_path)


def _confident_intervals(values: np.ndarray, threshold: float, min_length: int) -> list[str]:
    """Find contiguous confident runs at least ``min_length`` residues long."""
    mask = values >= threshold
    edges = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0] - 1
    lengths = ends - starts + 1
    return [
        f"{int(start) + 1}-{int(end) + 1}"
        for start, end, length in zip(starts, ends, lengths, strict=True)
        if length >= min_length
    ]


def _complement(intervals: list[str], protein_length: int) -> list[str]:
    """Return the intervals not covered by ``intervals``."""
    result: list[str] = []
    cursor = 1
    for interval in intervals:
        start, end = (int(value) for value in interval.split("-"))
        if cursor < start:
            result.append(f"{cursor}-{start - 1}")
        cursor = end + 1
    if cursor <= protein_length:
        result.append(f"{cursor}-{protein_length}")
    return result


def _confident_records(plddt: pd.DataFrame, parameters: PipelineParameters) -> pd.DataFrame:
    """Build confident and non-confident segment records for every protein."""
    records = []
    for uniprot_id, values_text in zip(plddt["uniprot_id"], plddt["pLDDT"], strict=True):
        values = np.asarray(values_text.split(","), dtype=float)
        confident = merge_intervals(
            _confident_intervals(
                values, parameters.plddt_threshold, parameters.confident_min_length
            ),
            gap=parameters.confident_merge_gap,
        )
        non_confident = _complement(confident, len(values))
        records.extend(
            {"uniprot_id": uniprot_id, "domain_residues": interval, "domain_type": "confident"}
            for interval in confident
        )
        records.extend(
            {"uniprot_id": uniprot_id, "domain_residues": interval, "domain_type": ""}
            for interval in non_confident
        )
    return pd.DataFrame.from_records(records)


def build_confident_table(
    plddt_path: Path,
    gene_path: Path,
    parameters: PipelineParameters,
    output_path: Path,
) -> None:
    """Build the 3D-AF-Surfer style confident-domain table."""
    plddt = read_tsv(plddt_path)
    genes = read_tsv(gene_path)
    require_columns(plddt, PLDDT_REQUIRED, plddt_path)
    require_columns(genes, GENE_REQUIRED, gene_path)

    records = _confident_records(plddt, parameters)
    records["_start"] = records["domain_residues"].str.split("-").str[0].astype(int)
    records = records.sort_values(["uniprot_id", "_start"], kind="stable").drop(columns="_start")
    records["domain_id"] = records.groupby("uniprot_id", sort=False).cumcount() + 1

    lookup = build_plddt_lookup(plddt)
    records = records.merge(plddt[["uniprot_id", "protein_length"]], on="uniprot_id", how="left")
    records = attach_stats(records, lookup, parameters.plddt_threshold)
    annotated = _join_gene_annotation(records, genes, fill_missing=False)
    write_tsv(annotated[CONFIDENT_COLUMNS], output_path)
