"""
AFdb structure inventory.
=========================

Scans the AlphaFold DB release once and writes the per-entry inventory every
later stage reads: modelled sequence, description, length, CIF filename, and
PAE availability. The scan is deliberately separated from the comparison so
the expensive CIF parsing is cached as its own Snakemake rule.

Input
-----
- AFdb CIF directory (``AF-*-model_v*.cif.gz``).
- AFdb PAE directory (``AF-*-predicted_aligned_error_v*.json``).

Output
------
- One row per CIF entry, sorted by accession, with PAE availability joined
  from the PAE directory. PAE files without a CIF entry are reported and
  ignored.

Usage
-----
    from src.afdb import join_entry_table, read_cif_entries, read_pae_index
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from enum import StrEnum
import gzip
from pathlib import Path
import re
from typing import TextIO

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from loguru import logger
import pandas as pd

# =============================================================================
# 2. CONSTANTS & ENUMS
# =============================================================================
CIF_SEQUENCE_KEY = "_entity_poly.pdbx_seq_one_letter_code_can"
CIF_DESCRIPTION_KEY = "_entity.pdbx_description"
CIF_NAME_PATTERN = re.compile(r"^AF-(?P<accession>.+)-F\d+-model_v\d+\.cif\.gz$")
PAE_NAME_PATTERN = re.compile(r"^AF-(?P<accession>.+)-F\d+-predicted_aligned_error_v\d+\.json$")
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


class AfdbCol(StrEnum):
    """Column names of the AFdb inventory tables."""

    ACCESSION = "afdb_accession"
    DESCRIPTION = "afdb_description"
    LENGTH = "afdb_length"
    SEQUENCE = "afdb_sequence"
    CIF_FILE = "cif_file"
    PAE_AVAILABLE = "pae_available"
    PAE_FILE = "pae_file"


CIF_ENTRY_COLUMNS = [
    AfdbCol.ACCESSION,
    AfdbCol.DESCRIPTION,
    AfdbCol.LENGTH,
    AfdbCol.SEQUENCE,
    AfdbCol.CIF_FILE,
]
PAE_INDEX_COLUMNS = [AfdbCol.ACCESSION, AfdbCol.PAE_FILE]
ENTRY_TABLE_COLUMNS = [*CIF_ENTRY_COLUMNS, AfdbCol.PAE_AVAILABLE, AfdbCol.PAE_FILE]

# =============================================================================
# 3. BOUNDARY LOADERS
# =============================================================================


def read_cif_header_sequence(handle: TextIO) -> tuple[str, str]:
    """Read the entity description and modelled sequence from an AFdb CIF header."""
    cif = MMCIF2Dict(handle)
    description = _first_cif_value(cif.get(CIF_DESCRIPTION_KEY, ""))
    sequence = _first_cif_value(cif.get(CIF_SEQUENCE_KEY, ""))
    return description, "".join(sequence.split())


def _first_cif_value(value: str | list[str]) -> str:
    """Return the first value of a CIF tag that may repeat across entities."""
    return value[0] if isinstance(value, list) else value


def parse_afdb_cif(path: Path) -> dict[str, str]:
    """Extract accession, description, and modelled sequence from one AFdb CIF."""
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
        AfdbCol.CIF_FILE: path.name,
    }


def read_cif_entries(afdb_dir: Path) -> pd.DataFrame:
    """Read every AFdb CIF header into one inventory row per entry."""
    paths = sorted(afdb_dir.glob("AF-*-F*-model_v*.cif.gz"))
    if not paths:
        raise FileNotFoundError(f"No AFdb CIF files found under {afdb_dir}")

    entries = pd.DataFrame.from_records([parse_afdb_cif(path) for path in paths])
    entries[AfdbCol.LENGTH] = entries[AfdbCol.SEQUENCE].str.len()

    unexpected = set("".join(entries[AfdbCol.SEQUENCE])) - STANDARD_AMINO_ACIDS
    if unexpected:
        logger.warning(f"Non-standard residues in AFdb sequences: {sorted(unexpected)}")

    entries = entries.sort_values(AfdbCol.ACCESSION, ignore_index=True)
    logger.info(f"Read {len(entries):,} AFdb CIF entries from {afdb_dir}")
    return entries.loc[:, CIF_ENTRY_COLUMNS]


def read_pae_index(pae_dir: Path) -> pd.DataFrame:
    """List the accession and filename of every AFdb PAE matrix."""
    paths = sorted(pae_dir.glob("AF-*-F*-predicted_aligned_error_v*.json"))
    if not paths:
        raise FileNotFoundError(f"No AFdb PAE files found under {pae_dir}")

    records: list[dict[str, str]] = []
    for path in paths:
        name_match = PAE_NAME_PATTERN.match(path.name)
        if name_match is None:
            raise ValueError(f"Unrecognised AFdb PAE filename: {path.name}")
        records.append(
            {AfdbCol.ACCESSION: name_match.group("accession"), AfdbCol.PAE_FILE: path.name}
        )

    index = pd.DataFrame.from_records(records).sort_values(AfdbCol.ACCESSION, ignore_index=True)
    logger.info(f"Indexed {len(index):,} AFdb PAE files from {pae_dir}")
    return index.loc[:, PAE_INDEX_COLUMNS]


# =============================================================================
# 4. FUNCTIONAL CORE
# =============================================================================
def join_entry_table(entries: pd.DataFrame, pae_index: pd.DataFrame) -> pd.DataFrame:
    """Join PAE availability onto the CIF inventory, reporting PAE-only accessions."""
    table = entries.merge(pae_index, on=AfdbCol.ACCESSION, how="left", validate="one_to_one")
    table[AfdbCol.PAE_AVAILABLE] = table[AfdbCol.PAE_FILE].notna()
    table[AfdbCol.PAE_FILE] = table[AfdbCol.PAE_FILE].fillna("")

    pae_only = set(pae_index[AfdbCol.ACCESSION]) - set(entries[AfdbCol.ACCESSION])
    if pae_only:
        logger.warning(f"PAE files without a CIF entry ({len(pae_only)}): {sorted(pae_only)}")

    logger.info(
        f"AFdb entry table: {len(table):,} entries, "
        f"{int(table[AfdbCol.PAE_AVAILABLE].sum()):,} with PAE"
    )
    return table.loc[:, ENTRY_TABLE_COLUMNS]
