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
    _match_prediction,
    _selection_reason,
    build_registry_stats,
)


def _merged() -> pd.DataFrame:
    """Build a minimal comparison frame with one sequence match and one ID fallback."""
    return pd.DataFrame(
        {
            "gene_systematic_id": ["geneA", "geneB"],
            "uniprot_id": ["A", "B"],
            "afdb_accession": ["A", "B"],
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
            "afdb_sequence": ["AAA", "CCC"],
            "pombase_accession": ["A", "B"],
        }
    )


def _merged_alternate() -> pd.DataFrame:
    """Build a frame where geneB is paired to a sequence-identical alternate entry."""
    return _merged().assign(
        afdb_accession=["A", "B_alt"],
        afdb_sequence=["AAA", "CCC"],
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
    matched = _match_prediction(merged, ["A", "B"], {"geneA": "AAA", "geneB": "TTT"})
    assert matched.to_dict("records") == [
        {"gene_systematic_id": "geneA", "structure_id": "A", "structure_sequence": "AAA"}
    ]


def test_match_af3_requires_sequence_equality() -> None:
    merged = _merged()
    matched = _match_prediction(merged, ["B"], {"B": "CCC"}, sequence_key="uniprot_id")
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


def test_build_registry_stats_counts_downstream_sources() -> None:
    registry = pd.DataFrame(
        {
            "in_analysis": [True, True, True, False],
            "structure_source": [VDULAB_SOURCE, AF3_SOURCE, VDULAB_SOURCE, "none"],
            "sequence_category": [
                "substitution_only",
                "not_compared",
                "not_compared",
                "not_compared",
            ],
            "category": ["shared", "shared", "pombase_only", "afdb_only"],
            "selection_reason": [
                VDULAB_SOURCE,
                AF3_SOURCE,
                VDULAB_SOURCE,
                "excluded_prediction_sequence_mismatch",
            ],
            "ted_ecod_coordinate_ok": [True, False, True, False],
        }
    )
    stats = build_registry_stats(registry).set_index("metric")["value"].to_dict()
    assert stats == {
        "in_analysis": 3,
        "nonidentical_shared": 1,
        "vdulab_total": 2,
        "vdulab_nonidentical": 1,
        "vdulab_long": 1,
        "af3_total": 1,
        "excluded_prediction_sequence_mismatch": 1,
        "ted_ecod_coordinate_ok": 2,
    }


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
