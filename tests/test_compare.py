"""Unit tests for the sequence-first AFdb/PomBase pairing."""

import importlib.util
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "compare_afdb_pombase.py"
SPEC = importlib.util.spec_from_file_location("compare_afdb_pombase", SCRIPT)
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)


def _genes() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_systematic_id": ["g1", "g2", "g3", "g4"],
            "gene_name": ["", "", "", ""],
            "gene_product": ["", "", "", ""],
            "gene_type": ["protein coding gene"] * 4,
            "uniprot_id": ["A", "B", None, "E"],
        }
    )


def _peptides() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_systematic_id": ["g1", "g2", "g3", "g4"],
            "pombase_sequence": ["AAA", "CCC", "DDD", "EEE"],
            "pombase_length": [3, 3, 3, 3],
            "internal_stop_truncated": [False, False, False, False],
        }
    )


def _entries() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "uniprot_id": ["A", "A_alt", "B", "C", "D"],
            "afdb_sequence": ["AAA", "AAA", "CCD", "CCC", "DDD"],
        }
    )


def test_sequence_match_prefers_own_accession() -> None:
    pairs = COMPARE.assign_structure_accessions(_genes(), _peptides(), _entries())
    row = pairs.set_index("gene_systematic_id").loc["g1"]
    assert row["match_method"] == "sequence"
    assert row["afdb_accession"] == "A"
    assert row["duplicate_afdb_accessions"] == "A_alt"


def test_sequence_match_uses_alternate_when_own_entry_differs() -> None:
    pairs = COMPARE.assign_structure_accessions(_genes(), _peptides(), _entries())
    row = pairs.set_index("gene_systematic_id").loc["g2"]
    assert row["match_method"] == "sequence"
    assert row["afdb_accession"] == "C"
    assert row["uniprot_id"] == "B"


def test_sequence_match_fills_missing_accession() -> None:
    pairs = COMPARE.assign_structure_accessions(_genes(), _peptides(), _entries())
    row = pairs.set_index("gene_systematic_id").loc["g3"]
    assert row["match_method"] == "sequence"
    assert row["afdb_accession"] == "D"
    assert row["uniprot_id"] == "D"


def test_unmatched_gene_stays_pombase_only() -> None:
    pairs = COMPARE.assign_structure_accessions(_genes(), _peptides(), _entries())
    row = pairs.set_index("gene_systematic_id").loc["g4"]
    assert row["match_method"] == "none"
    assert pd.isna(row["afdb_accession"])
    assert row["uniprot_id"] == "E"
