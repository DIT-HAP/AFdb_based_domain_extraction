"""
PAE matrix clustering.
======================

Parses AlphaFold predicted-aligned-error JSON files and partitions residues
into pseudo-rigid communities with the Leiden algorithm (python-igraph).

The random number generator is seeded explicitly so repeated runs produce the
same segmentation. This replaces the legacy ``fix_random_seed`` duplicate.

Input
-----
- A directory of ``*predicted_aligned_error*.json`` files.

Output
------
- One ``.domain`` file per input: comma-separated residue indices per cluster.

Usage
-----
    python scripts/afdb_domains.py pae-to-domains --pae-dir <dir> --output-dir <dir>
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import random

import igraph
import numpy as np
import pandas as pd

from src.files import AppConfig, PipelineParameters, read_tsv, require_columns, write_tsv
from src.registry import (
    AF3_SOURCE,
    AFDB_ALTERNATE,
    AFDB_PRIMARY,
    AFDB_SEQUENCE_MATCH,
    VDULAB_SOURCE,
)

# =============================================================================
# 2. CONSTANTS
# =============================================================================
DOMAIN_SUFFIX = ".domain"
VERSION_SEPARATOR = "-F1-"
REGISTRY_REQUIRED = ["gene_systematic_id", "in_analysis", "pae_path", "structure_source"]
SOURCE_TAGS = {
    AFDB_PRIMARY: "AFDBv6",
    AFDB_SEQUENCE_MATCH: "AFDBv6",
    AFDB_ALTERNATE: "AFDBv6",
    VDULAB_SOURCE: "vDulab",
    AF3_SOURCE: "vAF3",
}


# =============================================================================
# 3. FUNCTIONAL CORE
# =============================================================================
def parse_pae(path: Path) -> np.ndarray:
    """Parse one AlphaFold PAE JSON file into a square matrix."""
    payload = json.loads(path.read_text())[0]
    if "predicted_aligned_error" in payload:
        return np.asarray(payload["predicted_aligned_error"], dtype=np.float64)
    if "residue1" in payload and "distance" in payload:
        size = int(max(payload["residue1"]))
        matrix = np.empty((size, size), dtype=np.float64)
        matrix.ravel()[:] = payload["distance"]
        return matrix
    raise ValueError(f"Invalid PAE JSON format: {path}")


def cluster_pae_matrix(
    matrix: np.ndarray,
    *,
    pae_power: float,
    pae_cutoff: float,
    resolution: float,
) -> list[list[int]]:
    """Partition a PAE matrix into residue communities with Leiden clustering."""
    weights = 1 / matrix**pae_power
    edges = np.argwhere(matrix < pae_cutoff)
    graph = igraph.Graph(n=weights.shape[0], edges=edges.tolist(), directed=False)
    graph.es["weight"] = weights[edges[:, 0], edges[:, 1]]

    random.seed(0)
    igraph.set_random_number_generator(random)
    partition = graph.community_leiden(
        weights="weight",
        resolution_parameter=resolution / 100,
        n_iterations=-1,
    )

    members: dict[int, list[int]] = defaultdict(list)
    for residue, cluster in enumerate(partition.membership):
        members[int(cluster)].append(residue)
    return sorted(members.values(), key=len, reverse=True)


def write_domain_file(clusters: list[list[int]], path: Path) -> None:
    """Write one padded comma-separated cluster file."""
    width = max(len(cluster) for cluster in clusters)
    rows = [
        ",".join([*(str(residue) for residue in cluster), *([""] * (width - len(cluster)))])
        for cluster in clusters
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + "\n")


def domain_file_to_runs(path: Path) -> list[str]:
    """Convert a cluster file into comma-separated 1-based residue runs."""
    runs: list[str] = []
    for line in path.read_text().splitlines():
        residues = np.asarray([int(value) for value in line.split(",") if value], dtype=int)
        if residues.size == 0:
            continue
        breaks = np.where(np.diff(residues) != 1)[0] + 1
        runs.append(
            ",".join(f"{part[0] + 1}-{part[-1] + 1}" for part in np.split(residues, breaks))
        )
    return runs


def parse_pae_id(name: str) -> str:
    """Extract the protein identifier from a PAE JSON file name."""
    return name.split(VERSION_SEPARATOR)[0].removeprefix("AF-")


def read_subset(path: Path) -> set[str]:
    """Read a one-identifier-per-line subset file."""
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def parse_domain_filename(name: str) -> tuple[str, str]:
    """Extract the UniProt accession and model version from a ``.domain`` name."""
    stem = name.removesuffix(DOMAIN_SUFFIX)
    uniprot_id = stem.split(VERSION_SEPARATOR)[0].removeprefix("AF-")
    version = stem.rsplit("_", 1)[1]
    return uniprot_id, version


def build_domain_table(domains_dir: Path, output_path: Path, id_column: str = "uniprot_id") -> None:
    """Collect every ``.domain`` file into the PAE domain table."""
    records: list[dict[str, object]] = []
    for path in sorted(domains_dir.glob(f"*{DOMAIN_SUFFIX}")):
        uniprot_id, version = parse_domain_filename(path.name)
        runs = domain_file_to_runs(path)
        if not runs:
            continue
        protein_length = max(int(run.split("-")[-1]) for run in runs)
        records.extend(
            {
                "uniprot_id": uniprot_id,
                "version": version,
                "domain_id": domain_id,
                "domain_residues": run,
                "protein_length": protein_length,
            }
            for domain_id, run in enumerate(runs, start=1)
        )

    table = pd.DataFrame.from_records(
        records,
        columns=["uniprot_id", "version", "domain_id", "domain_residues", "protein_length"],
    )
    table = table.sort_values(["uniprot_id", "version", "domain_id"], ignore_index=True)
    write_tsv(table.rename(columns={"uniprot_id": id_column}), output_path)


def process_pae_file(
    pae_path: Path,
    output_dir: Path,
    parameters: PipelineParameters,
    output_name: str | None = None,
) -> None:
    """Cluster one PAE file and write its ``.domain`` output."""
    matrix = parse_pae(pae_path)
    clusters = cluster_pae_matrix(
        matrix,
        pae_power=parameters.pae_power,
        pae_cutoff=parameters.pae_cutoff,
        resolution=parameters.resolution,
    )
    output_path = output_dir / (output_name or f"{pae_path.stem}{DOMAIN_SUFFIX}")
    write_domain_file(clusters, output_path)


def run_pae_to_domains(
    pae_dir: Path,
    output_dir: Path,
    parameters: PipelineParameters,
    jobs: int,
    subset_path: Path | None = None,
) -> None:
    """Cluster every PAE JSON file in a directory, optionally in parallel."""
    pae_files = sorted(pae_dir.glob("*.json"))
    if subset_path is not None:
        subset = read_subset(subset_path)
        pae_files = [path for path in pae_files if parse_pae_id(path.name) in subset]
    if not pae_files:
        raise FileNotFoundError(f"No PAE JSON files found in {pae_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    if jobs <= 1:
        for path in pae_files:
            process_pae_file(path, output_dir, parameters)
        return
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        futures = [
            pool.submit(process_pae_file, path, output_dir, parameters) for path in pae_files
        ]
        for future in futures:
            future.result()


def run_pae_to_domains_from_registry(
    registry_path: Path,
    output_dir: Path,
    jobs: int,
    app_config: AppConfig,
) -> None:
    """Cluster the PAE matrix selected for every registry gene."""
    registry = read_tsv(registry_path)
    require_columns(registry, REGISTRY_REQUIRED, registry_path)
    selected = registry.loc[registry["in_analysis"]]
    output_names = (
        "AF-"
        + selected["gene_systematic_id"].astype(str)
        + "-F1-model_"
        + selected["structure_source"].map(SOURCE_TAGS)
        + DOMAIN_SUFFIX
    )
    tasks = list(zip(selected["pae_path"].map(app_config.resolve), output_names, strict=True))
    if not tasks:
        raise FileNotFoundError(f"No registry rows to cluster in {registry_path}")
    output_dir.mkdir(parents=True, exist_ok=True)

    if jobs <= 1:
        for pae_path, output_name in tasks:
            process_pae_file(pae_path, output_dir, app_config.parameters, output_name)
        return
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        futures = [
            pool.submit(process_pae_file, pae_path, output_dir, app_config.parameters, output_name)
            for pae_path, output_name in tasks
        ]
        for future in futures:
            future.result()
