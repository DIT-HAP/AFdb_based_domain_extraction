"""
Genomic projection of protein domains.
======================================

Maps protein-level domains onto the representative coding transcript and then
onto individual nucleotides and merged genomic intervals.

The nucleotide table is 0-based half-open (BED convention); the protein domain
table remains 1-based inclusive. Conversion happens through the accumulated
CDS length, exactly as in the legacy pipeline.

Input
-----
- PomBase GFF3, PeptideStats, and the annotated protein domain table.

Output
------
- Nucleotide bed, nucleotide-to-domain table, merged domain table, and BED.
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from pathlib import Path
import re

import numpy as np
import pandas as pd

from src.files import read_tsv, require_columns, write_tsv
from src.plddt import merge_intervals

# =============================================================================
# 2. CONSTANTS
# =============================================================================
GFF_COLUMNS = ["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attribute"]
SYSTEMATIC_ID_PATTERN = re.compile(r"ID=(\S+?)(?:$|(?:|\.\d(?::\S+|));)")
TRANSCRIPT_PATTERN = re.compile(r"Parent=(\S+?)$")
NUCLEOTIDE_REQUIRED = [
    "#Chr",
    "Start",
    "End",
    "Systematic ID",
    "Type",
    "Strand",
    "Length",
    "Transcript",
    "Start_Region",
    "End_Region",
    "Length_Region",
    "Accumulated_CDS_bases",
]
DOMAIN_REQUIRED = ["Systematic ID", "domain_id", "domain_residues", "domain_type"]
MAPPED_REQUIRED = [
    "#Chr",
    "Start",
    "End",
    "Systematic ID",
    "Type",
    "Strand",
    "Length",
    "Transcript",
]
GROUP_KEYS = [
    "#Chr",
    "Systematic ID",
    "Strand",
    "Transcript",
    "Start_Region",
    "End_Region",
    "Length_Region",
    "domain_id",
    "domain_residues",
    "domain_type",
]
BED_COLUMNS = [
    "#Chr",
    "Start_Region",
    "End_Region",
    "domain_residues",
    "score",
    "Strand",
    "Start",
    "End",
    "itemRgb",
]


# =============================================================================
# 3. FUNCTIONAL CORE
# =============================================================================
def parse_gff(gff_path: Path) -> pd.DataFrame:
    """Parse a GFF3 file and extract systematic ids and transcripts."""
    gff = pd.read_csv(
        gff_path, sep="\t", comment="#", header=None, names=GFF_COLUMNS, low_memory=False
    )
    gff["Systematic ID"] = gff["Attribute"].str.extract(SYSTEMATIC_ID_PATTERN, expand=False)
    gff["Transcript"] = gff["Attribute"].str.extract(TRANSCRIPT_PATTERN, expand=False)
    return gff


def representative_transcripts(gff: pd.DataFrame, peptide_stats: pd.DataFrame) -> pd.DataFrame:
    """Keep transcripts whose CDS length matches the PomBase peptide length."""
    coding_genes = gff.loc[gff["Feature"] == "mRNA", "Systematic ID"].unique()
    coding_transcripts = gff.loc[
        (gff["Feature"] == "CDS") & gff["Systematic ID"].isin(coding_genes),
        "Transcript",
    ].unique()
    coding_cds = gff.loc[
        (gff["Feature"] == "CDS") & gff["Transcript"].isin(coding_transcripts)
    ].copy()
    coding_cds["_span"] = coding_cds["End"] - coding_cds["Start"] + 1
    spans = (
        coding_cds.groupby(["Systematic ID", "Transcript"], sort=False)["_span"].sum().reset_index()
    )
    spans["Protein_length"] = (spans["_span"] - 3) / 3

    merged = spans.merge(
        peptide_stats[["Systematic_ID", "Residues"]],
        left_on="Systematic ID",
        right_on="Systematic_ID",
        how="outer",
    )
    merged["Delta"] = merged["Residues"] - merged["Protein_length"]
    return merged.loc[merged["Delta"] == 0].copy()


def cds_intron_table(gff: pd.DataFrame, peptide_stats: pd.DataFrame) -> pd.DataFrame:
    """Build the CDS/intron table of representative coding transcripts."""
    representatives = representative_transcripts(gff, peptide_stats)
    transcripts = representatives["Transcript"].dropna().unique()
    regions = gff.loc[
        gff["Transcript"].isin(transcripts) & gff["Feature"].isin(["CDS", "intron"]),
        ["Chr", "Start", "End", "Strand", "Feature", "Systematic ID", "Transcript"],
    ].rename(columns={"Feature": "Type"})
    regions["Start"] = regions["Start"] - 1
    regions["Length"] = regions["End"] - regions["Start"]

    coding = gff.loc[gff["Transcript"].isin(transcripts) & (gff["Feature"] == "CDS")]
    gene_region = (
        coding.groupby(["Systematic ID"], sort=False)
        .agg(
            Start_Region=("Start", lambda values: int(values.min()) - 1), End_Region=("End", "max")
        )
        .reset_index()
    )
    gene_region["Length_Region"] = gene_region["End_Region"] - gene_region["Start_Region"]

    merged = regions.merge(gene_region, on="Systematic ID", how="left")
    kept = merged.loc[
        (merged["Start"] >= merged["Start_Region"]) & (merged["End"] <= merged["End_Region"])
    ]
    return (
        kept[
            [
                "Chr",
                "Start",
                "End",
                "Systematic ID",
                "Type",
                "Strand",
                "Length",
                "Transcript",
                "Start_Region",
                "End_Region",
                "Length_Region",
            ]
        ]
        .rename(columns={"Chr": "#Chr"})
        .sort_values(["#Chr", "Start"], kind="stable", ignore_index=True)
    )


def expand_to_nucleotides(regions: pd.DataFrame) -> pd.DataFrame:
    """Expand half-open intervals to one row per nucleotide, fully vectorized."""
    lengths = regions["Length"].to_numpy(dtype=np.int64)
    total = int(lengths.sum())
    row_index = np.repeat(np.arange(len(regions)), lengths)
    offsets = np.arange(total, dtype=np.int64) - np.repeat(np.cumsum(lengths) - lengths, lengths)
    expanded = regions.iloc[row_index].reset_index(drop=True)
    expanded["Start"] = regions["Start"].to_numpy(dtype=np.int64)[row_index] + offsets
    expanded["End"] = expanded["Start"] + 1
    expanded["Length"] = 1
    return expanded


def accumulate_cds_bases(nucleotides: pd.DataFrame) -> pd.DataFrame:
    """Accumulate CDS bases in transcript order, strand-aware."""
    result = nucleotides.copy()
    result["_order"] = np.where(result["Strand"] == "-", -result["Start"], result["Start"])
    result = result.sort_values(["Systematic ID", "_order"], kind="stable")
    increment = result["Length"].where(result["Type"] == "CDS", 0)
    cumulative = increment.groupby(result["Systematic ID"], sort=False).cumsum()
    result["Accumulated_CDS_bases"] = cumulative - increment
    result = result.drop(columns="_order")
    return result.sort_values(["Systematic ID", "Start"], kind="stable", ignore_index=True)


# =============================================================================
# 4. IMPERATIVE SHELL
# =============================================================================
def build_nucleotide_bed(gff_path: Path, peptide_stats_path: Path, output_path: Path) -> None:
    """Create the per-nucleotide table for every coding region."""
    gff = parse_gff(gff_path)
    peptide_stats = read_tsv(peptide_stats_path)
    require_columns(peptide_stats, ["Systematic_ID", "Residues"], peptide_stats_path)
    regions = cds_intron_table(gff, peptide_stats)
    nucleotides = expand_to_nucleotides(regions)
    write_tsv(accumulate_cds_bases(nucleotides), output_path)


def map_domains_to_nucleotides(
    nucleotide_path: Path,
    annotated_path: Path,
    output_path: Path,
) -> None:
    """Assign each nucleotide the first protein domain covering its residue."""
    nucleotides = read_tsv(nucleotide_path)
    domains = read_tsv(annotated_path)
    require_columns(nucleotides, NUCLEOTIDE_REQUIRED, nucleotide_path)
    require_columns(domains, DOMAIN_REQUIRED, annotated_path)

    nucleotides["Residue_affected"] = nucleotides["Accumulated_CDS_bases"] // 3 + 1
    nucleotides["Residue_frame"] = nucleotides["Accumulated_CDS_bases"] % 3
    nucleotides = nucleotides.reset_index().rename(columns={"index": "_row"})

    exploded = domains.assign(_interval=domains["domain_residues"].str.split(",")).explode(
        "_interval"
    )
    bounds = exploded["_interval"].str.split("-", expand=True)
    exploded["domain_start"] = bounds[0].astype(int)
    exploded["domain_end"] = bounds[1].astype(int)

    candidates = nucleotides.merge(
        exploded[
            [
                "Systematic ID",
                "domain_id",
                "domain_residues",
                "domain_type",
                "domain_start",
                "domain_end",
            ]
        ],
        on="Systematic ID",
        how="inner",
    )
    in_domain = (candidates["Residue_affected"] >= candidates["domain_start"]) & (
        candidates["Residue_affected"] <= candidates["domain_end"]
    )
    matched = (
        candidates.loc[in_domain]
        .sort_values(["_row", "domain_id"], kind="stable")
        .drop_duplicates("_row", keep="first")
    )
    result = nucleotides.merge(
        matched[["_row", "domain_id", "domain_residues", "domain_type"]],
        on="_row",
        how="left",
    ).drop(columns="_row")
    write_tsv(result, output_path)


def organize_domains(input_path: Path, output_path: Path) -> None:
    """Merge nucleotide-level domain assignments into genomic intervals."""
    frame = read_tsv(input_path)
    require_columns(
        frame, [*MAPPED_REQUIRED, "domain_id", "domain_residues", "domain_type"], input_path
    )
    frame["Start_End"] = frame["Start"].astype(str) + "-" + frame["End"].astype(str)

    grouped = (
        frame.groupby(GROUP_KEYS, sort=True, dropna=True)["Start_End"].apply(list).reset_index()
    )
    grouped["Merged_Start_End"] = [
        ",".join(merge_intervals(values)) for values in grouped["Start_End"]
    ]
    grouped = grouped.drop(columns="Start_End")

    exploded = grouped.assign(_interval=grouped["Merged_Start_End"].str.split(",")).explode(
        "_interval"
    )
    bounds = exploded["_interval"].str.split("-", expand=True).astype(int)
    exploded["_start"] = bounds[0]
    exploded["_end"] = bounds[1]
    exploded["_length"] = exploded["_end"] - exploded["_start"]
    is_minus = exploded["Strand"] == "-"
    region_start = np.where(is_minus, exploded["End_Region"], exploded["Start_Region"])
    ordered_start = np.where(is_minus, exploded["_end"], exploded["_start"])
    ordered_end = np.where(is_minus, exploded["_start"], exploded["_end"])
    exploded["_ratio_start"] = (
        np.abs(ordered_start - region_start) / exploded["Length_Region"]
    ).round(2)
    exploded["_ratio_end"] = (np.abs(ordered_end - region_start) / exploded["Length_Region"]).round(
        2
    )
    exploded["_ratio"] = (
        exploded["_ratio_start"].astype(str) + "-" + exploded["_ratio_end"].astype(str)
    )

    merged = (
        exploded.groupby(GROUP_KEYS, sort=True)[["_length", "_ratio"]]
        .agg(
            Merged_Length=("_length", "sum"),
            domain_fraction=("_ratio", ",".join),
        )
        .reset_index()
    )
    result = grouped.merge(merged, on=GROUP_KEYS, how="left")
    write_tsv(result, output_path)


def build_genome_bed(input_path: Path, output_path: Path) -> None:
    """Expand merged genomic intervals into a BED9-style table."""
    frame = read_tsv(input_path)
    require_columns(
        frame,
        ["#Chr", "Start_Region", "End_Region", "domain_residues", "Strand", "Merged_Start_End"],
        input_path,
    )
    expanded = frame.assign(_interval=frame["Merged_Start_End"].str.split(",")).explode("_interval")
    bounds = expanded["_interval"].str.split("-", expand=True).astype(int)
    expanded["Start"] = bounds[0]
    expanded["End"] = bounds[1]
    expanded["score"] = "."
    expanded["itemRgb"] = "0,0,0"
    bed = expanded[BED_COLUMNS].sort_values(["#Chr", "Start"], kind="stable", ignore_index=True)
    write_tsv(bed, output_path)
