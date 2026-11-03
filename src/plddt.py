"""
pLDDT extraction and interval arithmetic.
==========================================

Reads per-residue confidence scores from AlphaFold PDB or mmCIF models and
provides the interval helpers used by the domain-merging stage.

Coordinates are 1-based and inclusive throughout the project.

Input
-----
- A directory of ``AF-<id>-F1-model_*.pdb.gz`` / ``.cif`` / ``.cif.gz`` files.

Output
------
- A TSV with ``uniprot_id``, ``protein_length``, and comma-separated ``pLDDT``.
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
import gzip
from pathlib import Path
from typing import TextIO

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
import numpy as np
import pandas as pd

from src.files import write_tsv

# =============================================================================
# 2. CONSTANTS
# =============================================================================
STRUCTURE_SUFFIXES = (".pdb.gz", ".cif.gz", ".cif", ".pdb")
STRUCTURE_PRIORITY = (".pdb.gz", ".cif", ".cif.gz")
IDENTIFIER_SEPARATOR = "-F1-"


# =============================================================================
# 3. FUNCTIONAL CORE
# =============================================================================
def extract_plddt(path: Path) -> list[float]:
    """Extract per-residue pLDDT values (CA B-factor) from a structure file."""
    if path.name.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            structure = _parse_structure(path, handle)
    else:
        with path.open("rt") as handle:
            structure = _parse_structure(path, handle)
    return [
        float(residue["CA"].bfactor) for residue in structure.get_residues() if residue.has_id("CA")
    ]


def _parse_structure(path: Path, handle: TextIO) -> Structure:
    """Parse a PDB or mmCIF handle with the matching Biopython parser."""
    stem = path.name.removesuffix(".gz")
    if stem.endswith(".pdb"):
        return PDBParser().get_structure(stem, handle)
    if stem.endswith(".cif"):
        return MMCIFParser().get_structure(stem, handle)
    raise ValueError(f"Unknown structure format: {path.name}")


def parse_structure_id(name: str) -> str:
    """Extract the protein identifier from an AlphaFold structure file name."""
    stem = name
    for suffix in STRUCTURE_SUFFIXES:
        if stem.endswith(suffix):
            stem = stem.removesuffix(suffix)
            break
    return stem.split(IDENTIFIER_SEPARATOR)[0].removeprefix("AF-")


def collect_structure_files(structure_dir: Path) -> dict[str, Path]:
    """Collect one structure file per protein, preferring PDB over mmCIF."""
    selected: dict[str, Path] = {}
    for suffix in STRUCTURE_PRIORITY:
        for path in sorted(structure_dir.glob(f"*{suffix}")):
            selected.setdefault(parse_structure_id(path.name), path)
    return selected


def build_plddt_table(
    structure_dir: Path,
    output_path: Path,
    subset_path: Path | None = None,
) -> None:
    """Build the per-protein pLDDT table for a structure directory."""
    structure_files = collect_structure_files(structure_dir)
    if subset_path is not None:
        subset = {line.strip() for line in subset_path.read_text().splitlines() if line.strip()}
        structure_files = {
            uniprot_id: path for uniprot_id, path in structure_files.items() if uniprot_id in subset
        }
    if not structure_files:
        raise FileNotFoundError(f"No structure files found in {structure_dir}")

    records = []
    for uniprot_id, path in structure_files.items():
        values = extract_plddt(path)
        records.append(
            {
                "uniprot_id": uniprot_id,
                "protein_length": len(values),
                "pLDDT": ",".join(str(value) for value in values),
            }
        )

    table = pd.DataFrame.from_records(records).sort_values("uniprot_id", ignore_index=True)
    write_tsv(table, output_path)


def interval_bounds(interval_string: str) -> tuple[np.ndarray, np.ndarray]:
    """Split a ``start-end[,start-end]`` string into start and end arrays."""
    parts = [part.split("-") for part in interval_string.split(",")]
    starts = np.asarray([int(part[0]) for part in parts], dtype=int)
    ends = np.asarray([int(part[1]) for part in parts], dtype=int)
    return starts, ends


def merge_intervals(intervals: list[str], gap: int = 1) -> list[str]:
    """Merge overlapping or nearby intervals, preserving ascending order."""
    bounds = sorted(
        (int(start), int(end)) for start, end in (item.split("-") for item in intervals)
    )
    merged: list[list[int]] = []
    for start, end in bounds:
        if not merged or start > merged[-1][1] + gap:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [f"{start}-{end}" for start, end in merged]


def check_overlap(interval: str, interval_list: list[str]) -> bool:
    """Return whether an interval overlaps any interval of a list."""
    start1, end1 = (int(value) for value in interval.split("-"))
    for other in interval_list:
        start2, end2 = (int(value) for value in other.split("-"))
        contains = start1 <= start2 <= end2 <= end1 or start2 <= start1 <= end1 <= end2
        bridges = start1 <= start2 < end1 <= end2 or start2 <= start1 < end2 <= end1
        if contains or bridges:
            return True
    return False


def domain_stats(
    interval_string: str,
    plddt_values: np.ndarray,
    threshold: float,
) -> tuple[int, float, float]:
    """Compute length, mean pLDDT, and fraction above threshold for a domain."""
    starts, ends = interval_bounds(interval_string)
    values = np.concatenate(
        [plddt_values[start - 1 : end] for start, end in zip(starts, ends, strict=True)]
    )
    length = int(values.size)
    mean = round(float(values.mean()), 2)
    fraction = round(float((values >= threshold).mean()), 2)
    return length, mean, fraction
