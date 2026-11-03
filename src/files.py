"""
Project I/O boundary.
=====================

Loads ``config/config.yaml``, resolves the ``cache:`` / ``legacy:`` / ``raw:``
/ ``project:`` path grammar, and provides fail-fast TSV access.

The path grammar keeps machine-specific roots in exactly one place, so moving
the project to another server means editing four values at the top of the YAML.

Input
-----
- ``config/config.yaml`` with path roots, resource paths, and parameters.

Output
------
- Immutable :class:`AppConfig` consumed by every pipeline stage.
"""

# =============================================================================
# 1. IMPORTS
# =============================================================================
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any

from loguru import logger
import pandas as pd
import yaml

# =============================================================================
# 2. CONFIGURATION OBJECTS
# =============================================================================


@dataclass(kw_only=True, frozen=True)
class ResourcePaths:
    """Resolved paths and release labels for external resources."""

    pombase_gene_table: Path
    pombase_peptide_stats: Path
    pombase_gff3: Path
    pombase_features_gene_table: Path
    pombase_features_peptide_stats: Path
    ted_dump: Path
    ted_segmentation_dump: Path
    ted_release: str
    ecod_dump: Path
    ecod_release: str


@dataclass(kw_only=True, frozen=True)
class PipelineParameters:
    """Immutable numerical parameters shared by the pipeline stages."""

    pae_power: float
    pae_cutoff: float
    resolution: float
    plddt_threshold: float
    merge_gap: int
    confident_min_length: int
    confident_merge_gap: int


def _resolve_spec(spec: str, roots: dict[str, Path]) -> Path:
    """Resolve a prefixed path specification against a root mapping."""
    prefix, separator, relative = spec.partition(":")
    if separator and prefix in roots:
        return (roots[prefix] / relative).resolve()
    return (roots["project"] / spec).resolve()


@dataclass(kw_only=True, frozen=True)
class AppConfig:
    """Immutable runtime configuration container."""

    cache_root: Path
    legacy_root: Path
    raw_root: Path
    project_root: Path
    resources: ResourcePaths
    parameters: PipelineParameters


# =============================================================================
# 3. IMPERATIVE SHELL (Boundary Loading & TSV Access)
# =============================================================================


def setup_logger(verbose: bool = False) -> None:
    """Configure the shared Loguru sink used by every command."""
    logger.remove()
    level = "DEBUG" if verbose else "INFO"
    logger.add(
        sys.stderr,
        level=level,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level:<8}</level> | {message}",
    )


def load_config(config_path: Path) -> AppConfig:
    """Load and validate the project configuration file."""
    raw = yaml.safe_load(config_path.read_text())
    paths = raw["paths"]
    resources = raw["resources"]
    parameters = raw["parameters"]
    roots = {
        "cache": Path(paths["cache_root"]),
        "legacy": Path(paths["legacy_root"]),
        "raw": Path(paths["raw_root"]),
        "project": Path(paths["project_root"]).resolve(),
    }
    pombase = resources["pombase"]
    return AppConfig(
        cache_root=roots["cache"],
        legacy_root=roots["legacy"],
        raw_root=roots["raw"],
        project_root=roots["project"],
        resources=ResourcePaths(
            pombase_gene_table=_resolve_spec(pombase["gene_table"], roots),
            pombase_peptide_stats=_resolve_spec(pombase["peptide_stats"], roots),
            pombase_gff3=_resolve_spec(pombase["gff3"], roots),
            pombase_features_gene_table=_resolve_spec(pombase["features_gene_table"], roots),
            pombase_features_peptide_stats=_resolve_spec(pombase["features_peptide_stats"], roots),
            ted_dump=_resolve_spec(resources["ted"]["filtered_dump"], roots),
            ted_segmentation_dump=_resolve_spec(resources["ted"]["segmentation_dump"], roots),
            ted_release=resources["ted"]["release"],
            ecod_dump=_resolve_spec(resources["ecod"]["dump"], roots),
            ecod_release=resources["ecod"]["release"],
        ),
        parameters=PipelineParameters(
            pae_power=float(parameters["pae_power"]),
            pae_cutoff=float(parameters["pae_cutoff"]),
            resolution=float(parameters["resolution"]),
            plddt_threshold=float(parameters["plddt_threshold"]),
            merge_gap=int(parameters["merge_gap"]),
            confident_min_length=int(parameters["confident_min_length"]),
            confident_merge_gap=int(parameters["confident_merge_gap"]),
        ),
    )


def read_tsv(path: Path, **kwargs: Any) -> pd.DataFrame:
    """Read a tab-separated file with the project defaults."""
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def write_tsv(frame: pd.DataFrame, path: Path, **kwargs: Any) -> None:
    """Write a tab-separated file, creating its parent directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, **kwargs)


def require_columns(frame: pd.DataFrame, required: list[str], source: Path) -> None:
    """Fail fast when a boundary input is missing required columns."""
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"Missing required columns in {source}: {missing}")
