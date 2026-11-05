"""
PomBase subsets of TED and ECOD.
================================

The upstream TED dump covers all eukaryotes and is hundreds of gigabytes; the
ECOD dump is the eukaryotic subset. Both builders stream the raw dump once and
keep only the proteins in the PomBase gene table, producing small resources
that can be committed.

Input
-----
- Raw TED dump (21 columns, no header) and ECOD domain dump (5 columns) from
  the reference root.
- PomBase gene table providing the UniProt identifiers to keep.

Output
------
- ``resources/ted/pombe_ted_segments.tsv``
- ``resources/ecod/pombe_ecod_segments.tsv``
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from pathlib import Path
import re

import pandas as pd

from src.files import load_gene_table, write_tsv

# =============================================================================
# 2. CONSTANTS
# =============================================================================
TED_COLUMNS = [
    "ted_id",
    "md5_domain",
    "consensus_level",
    "chopping",
    "nres_domain",
    "num_segments",
    "plddt",
    "num_helix_strand_turn",
    "num_helix",
    "num_strand",
    "num_helix_strand",
    "num_turn",
    "proteome_id",
    "cath_label",
    "cath_assignment_level",
    "cath_assignment_method",
    "packing_density",
    "norm_rg",
    "tax_common_name",
    "tax_scientific_name",
    "tax_lineage",
]
TED_OUTPUT_COLUMNS = [
    "ted_id",
    "uniprot_id",
    "consensus_level",
    "chopping",
    "nres_domain",
    "num_segments",
    "plddt",
    "packing_density",
    "norm_rg",
    "cath_label",
]
ECOD_COLUMNS = ["Protein", "Column2", "Range", "Source", "ID"]
ECOD_OUTPUT_COLUMNS = ["Protein", "Column2", "Range", "Source", "ID", "uniprot_id"]
TED_ID_PATTERN = re.compile(r"AF-([A-Za-z0-9]+)-F")


# =============================================================================
# 3. FUNCTIONAL CORE
# =============================================================================
def load_pombase_uniprot_ids(gene_path: Path) -> set[str]:
    """Load the non-empty UniProt identifiers of PomBase protein-coding genes."""
    genes = load_gene_table(gene_path)
    return set(genes["uniprot_id"].dropna().astype(str))


def stream_ted_matches(ted_path: Path, keep_ids: set[str]) -> list[list[str]]:
    """Stream the TED dump and keep rows belonging to PomBase proteins."""
    kept: list[list[str]] = []
    pattern = TED_ID_PATTERN
    with ted_path.open("rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            match = pattern.match(fields[0])
            if match and match.group(1) in keep_ids:
                kept.append(fields)
    return kept


def build_ted_table(ted_path: Path, gene_path: Path, output_path: Path) -> None:
    """Build the PomBase subset of the TED domain dump."""
    keep_ids = load_pombase_uniprot_ids(gene_path)
    rows = stream_ted_matches(ted_path, keep_ids)
    frame = pd.DataFrame.from_records(rows, columns=TED_COLUMNS)
    frame["uniprot_id"] = frame["ted_id"].str.extract(TED_ID_PATTERN, expand=False)
    write_tsv(frame[TED_OUTPUT_COLUMNS], output_path)


def build_ecod_table(ecod_path: Path, gene_path: Path, output_path: Path) -> None:
    """Build the PomBase subset of the ECOD domain dump."""
    keep_ids = load_pombase_uniprot_ids(gene_path)
    ecod = pd.read_csv(ecod_path, sep=r"\s+", header=None, names=ECOD_COLUMNS, low_memory=False)
    ecod["uniprot_id"] = ecod["Protein"].str.removesuffix("_F1")
    filtered = ecod.loc[ecod["uniprot_id"].isin(keep_ids)]
    write_tsv(filtered[ECOD_OUTPUT_COLUMNS], output_path)
