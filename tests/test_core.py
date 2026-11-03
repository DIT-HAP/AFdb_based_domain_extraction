"""Unit tests for the functional core."""

from pathlib import Path

import numpy as np
import pandas as pd
from src.dna import organize_domains
from src.pae import domain_file_to_runs, parse_domain_filename
from src.plddt import check_overlap, domain_stats, merge_intervals, parse_structure_id
from src.segments import assign_type
from src.unify import expand_intervals


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
    name = "AF-P12345-F1-predicted_aligned_error_v4.domain"
    assert parse_domain_filename(name) == ("P12345", "v4")


def test_domain_file_to_runs(tmp_path: Path) -> None:
    path = tmp_path / "AF-P12345-F1-predicted_aligned_error_v4.domain"
    path.write_text("0,1,2,5,6,\n")
    assert domain_file_to_runs(path) == ["1-3,6-7"]


def test_parse_structure_id() -> None:
    assert parse_structure_id("AF-P12345-F1-model_v4.pdb.gz") == "P12345"
    assert parse_structure_id("AF-SPAC23G3.02c_1-702-F1-model_vAF3.cif") == "SPAC23G3.02c_1-702"


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
