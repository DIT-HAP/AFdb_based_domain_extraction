"""Unit tests for the functional core."""

from pathlib import Path

import numpy as np
import pandas as pd
from src.dna import organize_domains
from src.files import PipelineParameters, load_config, write_tsv
from src.pae import domain_file_to_runs, parse_domain_filename
from src.plddt import check_overlap, domain_stats, merge_intervals
from src.segments import (
    assign_type,
    build_confident_table_from_registry,
    build_reformated_frame,
)
from src.unify import attach_pombase_annotation, expand_intervals

PARAMETERS = PipelineParameters(
    pae_power=1.0,
    pae_cutoff=5.0,
    resolution=0.5,
    plddt_threshold=70.0,
    merge_gap=11,
    confident_min_length=50,
    confident_merge_gap=6,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "config" / "config.yaml"


def test_load_config_reads_datasets() -> None:
    config = load_config(CONFIG_PATH)
    assert config.datasets["LocalAF2_PomBaseDiffAFdb"].version == (
        "LocalAF2_PomBaseDiffAFdb_20240305"
    )
    assert config.datasets["AF3_LongProt"].version == "AF3_LongProt_20241101"
    assert config.datasets["AFDB_v4"].structure.startswith("reference:UP000002485")


def test_merge_intervals_adjacent() -> None:
    assert merge_intervals(["1-10", "11-20"], gap=1) == ["1-20"]
    assert merge_intervals(["1-10", "11-20"], gap=0) == ["1-10", "11-20"]


def test_merge_intervals_gap() -> None:
    assert merge_intervals(["1-10", "21-30"], gap=11) == ["1-30"]
    assert merge_intervals(["1-10", "21-30"], gap=10) == ["1-10", "21-30"]


def test_check_overlap() -> None:
    assert check_overlap("1-10", ["5-15"]) is True
    assert check_overlap("1-10", ["10-20"]) is False
    assert check_overlap("1-10", ["11-20"]) is False


def test_domain_stats() -> None:
    values = np.asarray([90.0] * 10 + [60.0] * 10)
    assert domain_stats("1-10", values, 70.0) == (10, 90.0, 1.0)
    assert domain_stats("1-20", values, 70.0) == (20, 75.0, 0.5)


def test_parse_domain_filename() -> None:
    legacy = "AF-P12345-F1-predicted_aligned_error_v4.domain"
    assert parse_domain_filename(legacy) == ("P12345", "v4")
    registry = "AF-P12345-F1-model_LocalAF2_PomBaseDiffAFdb.domain"
    assert parse_domain_filename(registry) == ("P12345", "LocalAF2_PomBaseDiffAFdb")


def test_domain_file_to_runs(tmp_path: Path) -> None:
    path = tmp_path / "AF-P12345-F1-predicted_aligned_error_v4.domain"
    path.write_text("0,1,2,5,6,\n")
    assert domain_file_to_runs(path) == ["1-3,6-7"]


def test_expand_intervals() -> None:
    frame = pd.DataFrame(
        {
            "uniprot_id": ["P1"],
            "segment_id": ["1"],
            "domain_residues": ["1-5,10-12"],
        }
    )
    expanded = expand_intervals(frame, "domain_residues", ",")
    assert expanded["start"].tolist() == [1, 10]
    assert expanded["end"].tolist() == [5, 12]
    assert expanded["length"].tolist() == [5, 3]
    assert expanded["segment_part"].tolist() == [1, 2]


def test_attach_pombase_annotation_fills_blank_gene_name() -> None:
    frame = pd.DataFrame(
        {
            "uniprot_id": ["P1", "P2"],
            "systematic_id": [pd.NA, pd.NA],
            "gene_name": [pd.NA, pd.NA],
            "protein_length": [100.0, np.nan],
        }
    )
    genes = pd.DataFrame(
        {
            "uniprot_id": ["P1", "P2"],
            "gene_systematic_id": ["SPAC1.01", "SPAC2.02"],
            "gene_name": ["abc", ""],
        }
    )
    stats = pd.DataFrame({"Systematic_ID": ["SPAC1.01", "SPAC2.02"], "Residues": [100, 200]})

    result = attach_pombase_annotation(frame, genes, stats)

    assert result["gene_name"].tolist() == ["abc", "SPAC2.02"]


def test_assign_type() -> None:
    frame = pd.DataFrame(
        {
            "domain_residues": ["1-50", "51-100", "1-100", "21-80"],
            "domain_length": [50, 50, 100, 60],
            "protein_length": [100, 100, 100, 100],
            "mean_pLDDT": [90.0, 60.0, 80.0, 60.0],
        }
    )
    result = assign_type(frame, 70.0)
    assert result["location"].tolist() == ["N-terminal", "C-terminal", "Whole", "middle"]
    assert result["domain_type"].tolist() == ["domain", "disordered", "domain", "disordered"]


def test_build_reformated_frame_uses_registry_gene_id() -> None:
    pae = pd.DataFrame(
        {
            "gene_systematic_id": ["g1", "g1"],
            "version": ["AF3_LongProt", "AF3_LongProt"],
            "protein_length": [6, 6],
            "domain_residues": ["1-3", "4-6"],
            "mean_pLDDT": [90.0, 50.0],
        }
    )
    plddt = pd.DataFrame(
        {
            "gene_systematic_id": ["g1"],
            "protein_length": [6],
            "pLDDT": ["90,90,90,50,50,50"],
        }
    )
    result = build_reformated_frame(pae, plddt, PARAMETERS, id_column="gene_systematic_id")
    assert result["gene_systematic_id"].tolist() == ["g1", "g1"]
    assert result["domain_residues"].tolist() == ["1-3", "4-6"]
    assert "uniprot_id" not in result.columns


def test_build_confident_table_from_registry_uses_gene_id(tmp_path: Path) -> None:
    plddt = pd.DataFrame(
        {
            "gene_systematic_id": ["g1"],
            "protein_length": [120],
            "pLDDT": [",".join(["90"] * 60 + ["50"] * 60)],
        }
    )
    registry = pd.DataFrame({"gene_systematic_id": ["g1"], "uniprot_id": ["P1"]})
    plddt_path = tmp_path / "plddt.tsv"
    registry_path = tmp_path / "registry.tsv"
    output_path = tmp_path / "confident.tsv"
    write_tsv(plddt, plddt_path)
    write_tsv(registry, registry_path)

    build_confident_table_from_registry(plddt_path, registry_path, PARAMETERS, output_path)

    result = pd.read_csv(output_path, sep="\t", keep_default_na=False)
    assert result["Systematic ID"].tolist() == ["g1", "g1"]
    assert result["uniprot_id"].tolist() == ["P1", "P1"]
    assert result["domain_residues"].tolist() == ["1-60", "61-120"]
    assert result["domain_type"].tolist() == ["confident", ""]
    assert result["protein_length"].tolist() == [120, 120]


def test_organize_domains(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {
            "#Chr": ["chr1", "chr1"],
            "Systematic ID": ["G1", "G1"],
            "Type": ["CDS", "CDS"],
            "Strand": ["+", "+"],
            "Transcript": ["G1.1", "G1.1"],
            "Start_Region": [0, 0],
            "End_Region": [10, 10],
            "Length_Region": [10, 10],
            "domain_id": [1, 1],
            "domain_residues": ["1-4", "1-4"],
            "domain_type": ["domain", "domain"],
            "Start": [0, 1],
            "End": [1, 2],
            "Length": [1, 1],
        }
    )
    input_path = tmp_path / "mapped.tsv"
    output_path = tmp_path / "merged.tsv"
    frame.to_csv(input_path, sep="\t", index=False)
    organize_domains(input_path, output_path)
    result = pd.read_csv(output_path, sep="\t")
    assert result["Merged_Start_End"].tolist() == ["0-2"]
    assert result["Merged_Length"].tolist() == [2]
    assert result["domain_fraction"].tolist() == ["0.0-0.2"]
