"""Unit tests for the gene structure registry selection rules."""

import pandas as pd
from src.files import load_gene_table
from src.registry import (
    AF3_SOURCE,
    AFDB_ALTERNATE,
    AFDB_PRIMARY,
    AFDB_SEQUENCE_MATCH,
    VDULAB_SOURCE,
    _attach_structure,
    _match_af3,
    _match_vdulab,
    _selection_reason,
)


def _merged() -> pd.DataFrame:
    """Build a minimal comparison frame with one sequence match and one ID fallback."""
    return pd.DataFrame(
        {
            "gene_systematic_id": ["geneA", "geneB"],
            "uniprot_id": ["A", "B"],
            "afdb_accession": ["A", "B"],
            "afdb_entry_id": ["A", "B"],
            "category": ["shared", "shared"],
            "match_method": ["uniprot_id", "uniprot_id"],
            "sequence_category": ["identical", "indel_only"],
            "pombase_sequence": ["AAA", "CCC"],
            "pombase_length": [3, 3],
            "afdb_length": [3, 3],
            "length_delta": [0, 0],
            "n_substitutions": [0, 0],
            "n_afdb_extra": [0, 0],
            "n_pombase_extra": [0, 0],
            "identity_pct": [100.0, 50.0],
            "duplicate_afdb_accessions": [None, None],
            "notes": [None, None],
            "gene_name": ["a", "b"],
            "gene_product": ["pa", "pb"],
            "gene_type": ["protein coding gene", "protein coding gene"],
            "afdb_description": ["da", "db"],
            "pombase_accession": ["A", "B"],
        }
    )


def _merged_alternate() -> pd.DataFrame:
    """Build a frame where geneB is paired to a sequence-identical alternate entry."""
    return _merged().assign(
        afdb_accession=["A", "B_alt"],
        afdb_entry_id=["A", "B_alt"],
        match_method=["sequence", "sequence"],
        sequence_category=["identical", "identical"],
        afdb_length=[3, 3],
        length_delta=[0, 0],
        n_substitutions=[0, 0],
        n_afdb_extra=[0, 0],
        n_pombase_extra=[0, 0],
        identity_pct=[100.0, 100.0],
    )


def test_match_vdulab_requires_sequence_equality() -> None:
    merged = _merged()
    matched = _match_vdulab(merged, ["A", "B"], {"geneA": "AAA", "geneB": "TTT"})
    assert matched.to_dict("records") == [
        {"gene_systematic_id": "geneA", "structure_id": "A", "structure_sequence": "AAA"}
    ]


def test_match_af3_requires_sequence_equality() -> None:
    merged = _merged()
    matched = _match_af3(merged, ["B"], {"B": "CCC"})
    assert matched.to_dict("records") == [
        {"gene_systematic_id": "geneB", "structure_id": "B", "structure_sequence": "CCC"}
    ]


def test_selection_reason_prefers_ordinary_and_alternate_sequence_matches() -> None:
    empty = pd.DataFrame({"gene_systematic_id": []})
    reasons = _selection_reason(_merged_alternate(), empty, empty)
    assert reasons.tolist() == [AFDB_PRIMARY, AFDB_ALTERNATE]


def test_attach_structure_marks_missing_accession_as_sequence_match() -> None:
    registry = _merged_alternate().assign(
        pombase_accession=[None, "B"],
        selection_reason=[AFDB_PRIMARY, AFDB_ALTERNATE],
    )
    empty = pd.DataFrame(columns=["gene_systematic_id", "structure_id", "structure_sequence"])
    result = _attach_structure(
        registry,
        empty,
        empty,
        entry_sequences={"A": "AAA", "B_alt": "CCC"},
        v4_sequences={},
    )
    source = result.set_index("gene_systematic_id")["structure_source"]
    assert source.loc["geneA"] == AFDB_SEQUENCE_MATCH


def test_selection_reason_falls_back_to_vdulab_and_af3() -> None:
    merged = _merged()
    empty = pd.DataFrame({"gene_systematic_id": []})
    vdulab = pd.DataFrame({"gene_systematic_id": ["geneB"]})
    reasons = _selection_reason(merged, vdulab, empty)
    assert reasons.tolist() == [AFDB_PRIMARY, VDULAB_SOURCE]

    af3 = pd.DataFrame({"gene_systematic_id": ["geneB"]})
    reasons = _selection_reason(merged, empty, af3)
    assert reasons.tolist() == [AFDB_PRIMARY, AF3_SOURCE]


def test_attach_structure_assigns_source_and_paths() -> None:
    registry = _merged_alternate().assign(
        selection_reason=[AFDB_PRIMARY, AFDB_ALTERNATE],
        relationship=[None, None],
    )
    vdulab = pd.DataFrame(columns=["gene_systematic_id", "structure_id", "structure_sequence"])
    af3 = pd.DataFrame(columns=["gene_systematic_id", "structure_id", "structure_sequence"])
    result = _attach_structure(
        registry,
        vdulab,
        af3,
        entry_sequences={"A": "AAA", "B_alt": "CCC"},
        v4_sequences={"A": "AAA", "B": "CCD"},
    )
    source = result.set_index("gene_systematic_id")["structure_source"]
    assert source.loc["geneA"] == AFDB_PRIMARY
    assert source.loc["geneB"] == AFDB_ALTERNATE
    assert result.set_index("gene_systematic_id")["ted_ecod_coordinate_ok"].to_dict() == {
        "geneA": True,
        "geneB": False,
    }
    assert result.set_index("gene_systematic_id")["modeled_length"].loc["geneB"] == 3


def test_load_gene_table_accepts_external_id(tmp_path) -> None:
    path = tmp_path / "genes.tsv"
    path.write_text(
        "gene_systematic_id\tgene_name\tgene_type\texternal_id\n"
        "SPAC1\ta\tprotein coding gene\tP1\n"
        "SPAC2\tb\tlncRNA gene\tURS1\n"
    )
    table = load_gene_table(path)
    assert table["uniprot_id"].tolist() == ["P1"]
    assert table["gene_systematic_id"].tolist() == ["SPAC1"]
