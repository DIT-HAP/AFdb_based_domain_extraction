"""Unit tests for the AFdb inventory stage."""

import gzip
from pathlib import Path

import pytest
from src.afdb import AfdbCol, join_entry_table, parse_afdb_cif, read_cif_entries, read_pae_index

CIF_TEMPLATE = """data_{accession}
#
_entity.pdbx_description         "{description}"
#
_entity_poly.pdbx_seq_one_letter_code_can
;{sequence}
;
#
"""


def _write_cif(directory: Path, accession: str, sequence: str, description: str) -> None:
    """Write one minimal gzipped AFdb CIF header."""
    path = directory / f"AF-{accession}-F1-model_v6.cif.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(
            CIF_TEMPLATE.format(accession=accession, description=description, sequence=sequence)
        )


def _write_pae(directory: Path, accession: str) -> None:
    """Write one empty PAE JSON placeholder with the AFdb filename pattern."""
    (directory / f"AF-{accession}-F1-predicted_aligned_error_v6.json").write_text("{}")


def test_read_cif_entries_sorts_and_computes_length(tmp_path: Path) -> None:
    _write_cif(tmp_path, "P2", "MSTN", "Second")
    _write_cif(tmp_path, "P1", "MSTNPQ", "First")

    entries = read_cif_entries(tmp_path)

    assert entries[AfdbCol.ACCESSION].tolist() == ["P1", "P2"]
    assert entries[AfdbCol.LENGTH].tolist() == [6, 4]
    assert entries[AfdbCol.DESCRIPTION].tolist() == ["First", "Second"]
    assert entries[AfdbCol.CIF_FILE].tolist() == [
        "AF-P1-F1-model_v6.cif.gz",
        "AF-P2-F1-model_v6.cif.gz",
    ]


def test_parse_afdb_cif_rejects_empty_sequence(tmp_path: Path) -> None:
    path = tmp_path / "AF-P3-F1-model_v6.cif.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("data_P3\n#\n")

    with pytest.raises(ValueError, match="Empty sequence"):
        parse_afdb_cif(path)


def test_read_pae_index_lists_accessions(tmp_path: Path) -> None:
    _write_pae(tmp_path, "P9")
    _write_pae(tmp_path, "P1")

    index = read_pae_index(tmp_path)

    assert index[AfdbCol.ACCESSION].tolist() == ["P1", "P9"]
    assert index[AfdbCol.PAE_FILE].tolist() == [
        "AF-P1-F1-predicted_aligned_error_v6.json",
        "AF-P9-F1-predicted_aligned_error_v6.json",
    ]


def test_join_entry_table_keeps_entries_without_pae(tmp_path: Path) -> None:
    _write_cif(tmp_path, "P1", "MSTN", "First")
    _write_cif(tmp_path, "P2", "MSTK", "Second")
    _write_pae(tmp_path, "P1")
    _write_pae(tmp_path, "P9")

    table = join_entry_table(read_cif_entries(tmp_path), read_pae_index(tmp_path))

    assert table[AfdbCol.ACCESSION].tolist() == ["P1", "P2"]
    assert table[AfdbCol.PAE_AVAILABLE].tolist() == [True, False]
    assert table[AfdbCol.PAE_FILE].tolist() == [
        "AF-P1-F1-predicted_aligned_error_v6.json",
        "",
    ]
