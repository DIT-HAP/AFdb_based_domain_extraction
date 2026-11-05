"""Unit tests for the UniProt refresh assembly."""

from typing import Any

import pandas as pd
from src.uniprot import (
    AFDB_MATCHES,
    NEITHER_MATCHES,
    POMBASE_MATCHES,
    UniprotCol,
    UniprotEntry,
    build_resolution,
    build_triangulation,
    parse_entry,
)


def _entry(primary: str, sequence: str) -> UniprotEntry:
    """Build one active UniProt entry for the assembly tests."""
    return UniprotEntry(
        primary_accession=primary,
        reviewed="reviewed",
        protein_name=f"Protein {primary}",
        length=len(sequence),
        sequence=sequence,
    )


def test_parse_entry_reads_recommended_name() -> None:
    payload: dict[str, Any] = {
        "primaryAccession": "P1",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Protein one"}}},
        "sequence": {"length": 4, "value": "MSTN"},
    }

    assert parse_entry(payload) == UniprotEntry(
        primary_accession="P1",
        reviewed="reviewed",
        protein_name="Protein one",
        length=4,
        sequence="MSTN",
    )


def test_build_triangulation_labels_sequence_source() -> None:
    pairs = pd.DataFrame(
        {
            UniprotCol.GENE_SYSTEMATIC_ID: ["g1", "g2", "g3"],
            UniprotCol.UNIPROT_ID: ["P1", "P2", "P3"],
            UniprotCol.AFDB_SEQUENCE: ["MSTN", "MSTN", "MSTN"],
            UniprotCol.POMBASE_SEQUENCE: ["MSTK", "MSTK", "MSTK"],
        }
    )
    entries = {
        "P1": _entry("P1", "MSTK"),
        "P2": _entry("P2", "MSTN"),
        "P3": _entry("P3", "MSTM"),
    }

    table = build_triangulation(pairs, entries)

    assert table[UniprotCol.RELATIONSHIP].tolist() == [
        POMBASE_MATCHES,
        AFDB_MATCHES,
        NEITHER_MATCHES,
    ]
    assert table[UniprotCol.UNIPROT_LENGTH].tolist() == [4, 4, 4]


def test_build_resolution_marks_secondary_accessions() -> None:
    queries = pd.DataFrame({UniprotCol.UNIPROT_ID: ["A1", "A2"]})
    entries = {"A1": _entry("P1", "MSTN"), "A2": _entry("A2", "MSTK")}

    table = build_resolution(queries, entries, {"P1"})

    assert table[UniprotCol.QUERY_ACCESSION].tolist() == ["A1", "A2"]
    assert table[UniprotCol.IS_SECONDARY].tolist() == [True, False]
    assert table[UniprotCol.PRIMARY_ACCESSION].tolist() == ["P1", "A2"]
    assert table[UniprotCol.PRIMARY_IN_POMBASE].tolist() == [True, False]
