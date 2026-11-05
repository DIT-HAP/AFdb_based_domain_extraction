"""
UniProt snapshot refresh.
=========================

Fetches the UniProt records behind the two pinned snapshots in
``resources/uniprot/``: the non-identical PomBase/AFdb pairs and the AFdb-only
accessions. The refresh is a manual, explicit operation; the default pipeline
reads the committed snapshots so results stay reproducible.

Input
-----
- ``05_nonidentical_pairs.tsv`` and ``06_afdb_only.tsv`` from the comparison.
- PomBase gene table for the ``primary_in_pombase`` column.

Output
------
- Minimal-column refresh candidates under ``results/02_uniprot/``.

Usage
-----
    refresh-uniprot-triangulation --nonidentical-pairs ... --output ...
    refresh-uniprot-resolution --afdb-only ... --gene-table ... --output ...
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from dataclasses import dataclass
from enum import StrEnum
import time
from typing import Any

import httpx
from loguru import logger
import numpy as np
import pandas as pd

# =============================================================================
# 2. CONSTANTS & ENUMS
# =============================================================================
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb"
INACTIVE_ENTRY_TYPE = "Inactive"
REVIEWED_PREFIX = "UniProtKB reviewed"
REVIEWED_LABEL = "reviewed"
UNREVIEWED_LABEL = "unreviewed"
REQUEST_TIMEOUT_SECONDS = 30.0
REQUEST_DELAY_SECONDS = 1.0
MAX_ATTEMPTS = 3
RETRY_STATUS_CODES = {429, 500, 502, 503, 504}

POMBASE_MATCHES = "pombase_matches_uniprot_afdb_differs"
AFDB_MATCHES = "afdb_matches_uniprot_pombase_differs"
NEITHER_MATCHES = "neither_matches_uniprot"


class UniprotCol(StrEnum):
    """Column names of the UniProt refresh tables."""

    QUERY_ACCESSION = "query_accession"
    IS_SECONDARY = "is_secondary"
    PRIMARY_ACCESSION = "primary_accession"
    PRIMARY_IN_POMBASE = "primary_in_pombase"
    REVIEWED = "reviewed"
    LENGTH = "length"
    PROTEIN_NAME = "protein_name"
    GENE_SYSTEMATIC_ID = "gene_systematic_id"
    UNIPROT_ID = "uniprot_id"
    AFDB_SEQUENCE = "afdb_sequence"
    POMBASE_SEQUENCE = "pombase_sequence"
    RELATIONSHIP = "relationship"
    UNIPROT_LENGTH = "uniprot_length"


RESOLUTION_COLUMNS = [
    UniprotCol.QUERY_ACCESSION,
    UniprotCol.IS_SECONDARY,
    UniprotCol.PRIMARY_ACCESSION,
    UniprotCol.PRIMARY_IN_POMBASE,
    UniprotCol.REVIEWED,
    UniprotCol.LENGTH,
    UniprotCol.PROTEIN_NAME,
]
TRIANGULATION_COLUMNS = [
    UniprotCol.GENE_SYSTEMATIC_ID,
    UniprotCol.UNIPROT_ID,
    UniprotCol.RELATIONSHIP,
    UniprotCol.UNIPROT_LENGTH,
]
TRIANGULATION_INPUT_COLUMNS = [
    UniprotCol.GENE_SYSTEMATIC_ID,
    UniprotCol.UNIPROT_ID,
    UniprotCol.AFDB_SEQUENCE,
    UniprotCol.POMBASE_SEQUENCE,
]
RESOLUTION_INPUT_COLUMNS = [UniprotCol.UNIPROT_ID]

# =============================================================================
# 3. CONFIGURATION
# =============================================================================


@dataclass(frozen=True, kw_only=True)
class UniprotEntry:
    """One active UniProtKB entry resolved for a query accession."""

    primary_accession: str
    reviewed: str
    protein_name: str
    length: int
    sequence: str


# =============================================================================
# 4. FUNCTIONAL CORE
# =============================================================================
def _protein_name(payload: dict[str, Any]) -> str:
    """Read the recommended or first submitted protein name from a UniProtKB record."""
    description = payload.get("proteinDescription", {})
    recommended = description.get("recommendedName", {}).get("fullName", {}).get("value")
    if recommended:
        return str(recommended)
    submissions = description.get("submissionNames", [])
    return str(submissions[0]["fullName"]["value"]) if submissions else ""


def _reviewed_label(entry_type: str) -> str:
    """Map the UniProtKB entry type to the TSV reviewed label."""
    return REVIEWED_LABEL if entry_type.startswith(REVIEWED_PREFIX) else UNREVIEWED_LABEL


def parse_entry(payload: dict[str, Any]) -> UniprotEntry:
    """Parse the fields needed from an active UniProtKB JSON entry."""
    sequence = payload["sequence"]
    return UniprotEntry(
        primary_accession=str(payload["primaryAccession"]),
        reviewed=_reviewed_label(str(payload["entryType"])),
        protein_name=_protein_name(payload),
        length=int(sequence["length"]),
        sequence=str(sequence["value"]),
    )


def build_resolution(
    queries: pd.DataFrame,
    entries: dict[str, UniprotEntry],
    pombase_accessions: set[str],
) -> pd.DataFrame:
    """Build the minimal resolution table from fetched entries."""
    records = [
        {
            UniprotCol.QUERY_ACCESSION: accession,
            UniprotCol.IS_SECONDARY: entries[accession].primary_accession != accession,
            UniprotCol.PRIMARY_ACCESSION: entries[accession].primary_accession,
            UniprotCol.PRIMARY_IN_POMBASE: entries[accession].primary_accession
            in pombase_accessions,
            UniprotCol.REVIEWED: entries[accession].reviewed,
            UniprotCol.LENGTH: entries[accession].length,
            UniprotCol.PROTEIN_NAME: entries[accession].protein_name,
        }
        for accession in queries[UniprotCol.UNIPROT_ID]
    ]
    return pd.DataFrame.from_records(records, columns=RESOLUTION_COLUMNS)


def build_triangulation(
    pairs: pd.DataFrame,
    entries: dict[str, UniprotEntry],
) -> pd.DataFrame:
    """Rebuild the triangulation table from current sequences and fetched entries."""
    uniprot_sequences = pairs[UniprotCol.UNIPROT_ID].map(
        {accession: entry.sequence for accession, entry in entries.items()}
    )
    pombase_matches = uniprot_sequences == pairs[UniprotCol.POMBASE_SEQUENCE]
    afdb_matches = uniprot_sequences == pairs[UniprotCol.AFDB_SEQUENCE]

    relationship = np.select(
        [pombase_matches & ~afdb_matches, afdb_matches & ~pombase_matches],
        [POMBASE_MATCHES, AFDB_MATCHES],
        default=NEITHER_MATCHES,
    )
    table = pairs.loc[:, [UniprotCol.GENE_SYSTEMATIC_ID, UniprotCol.UNIPROT_ID]].assign(
        **{
            UniprotCol.RELATIONSHIP: relationship,
            UniprotCol.UNIPROT_LENGTH: uniprot_sequences.str.len().astype("Int64"),
        }
    )
    return table.loc[:, TRIANGULATION_COLUMNS]


# =============================================================================
# 5. IMPERATIVE SHELL (HTTP Boundary)
# =============================================================================
def _get_json(client: httpx.Client, accession: str) -> dict[str, Any]:
    """Fetch one UniProtKB JSON record with retry and exponential backoff."""
    url = f"{UNIPROT_ENTRY_URL}/{accession}.json"
    for attempt in range(1, MAX_ATTEMPTS + 1):
        response = client.get(url)
        if response.status_code == httpx.codes.OK:
            return response.json()
        if response.status_code not in RETRY_STATUS_CODES or attempt == MAX_ATTEMPTS:
            response.raise_for_status()
        logger.warning(f"UniProt {accession}: HTTP {response.status_code}, retry {attempt}")
        time.sleep(REQUEST_DELAY_SECONDS * 2**attempt)
    raise RuntimeError(f"UniProt request failed for {accession}")


def fetch_entry(client: httpx.Client, accession: str) -> UniprotEntry:
    """Resolve one query accession to its active UniProtKB entry."""
    payload = _get_json(client, accession)
    if payload.get("entryType") == INACTIVE_ENTRY_TYPE:
        targets = payload.get("inactiveReason", {}).get("mergeDemergeTo", [])
        if not targets:
            raise ValueError(f"UniProt accession {accession} is inactive without a merge target")
        payload = _get_json(client, str(targets[0]))
    return parse_entry(payload)


def fetch_entries(accessions: list[str]) -> dict[str, UniprotEntry]:
    """Fetch every query accession, resolving merged secondary accessions."""
    entries: dict[str, UniprotEntry] = {}
    with httpx.Client(timeout=REQUEST_TIMEOUT_SECONDS, follow_redirects=True) as client:
        for accession in accessions:
            entries[accession] = fetch_entry(client, accession)
            time.sleep(REQUEST_DELAY_SECONDS)
    logger.info(f"Fetched {len(entries):,} UniProt entries")
    return entries
