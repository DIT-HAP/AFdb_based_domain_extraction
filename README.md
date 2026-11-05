# AFdb-based domain extraction

Call protein domains from AlphaFold PAE matrices, merge them with per-residue
pLDDT confidence, and emit one unified segment table per dataset alongside the
PomBase subsets of TED and ECOD.

The v2 refactor keeps the original numerical behaviour where it was defined and
replaces the ad-hoc shell orchestration with a Snakemake workflow.

## Deliverables

- `results/02_gene_structure_registry.tsv` — one row per PomBase coding
  gene (plus the AFdb-only entries), recording the selected structure, PAE
  matrix, modelled sequence, coordinate basis, and analysis flag.
- `results/09_segments.tsv` — the canonical unified segment table with
  three sources: `PAE`, `TED`, `ECOD`.
- `results/01_comparison/` — the detailed AFdb-vs-PomBase tables, report,
  and Excel workbook.

`docs/results_index.md` numbers every output file, its row count, the producing
rule (`R#`), and the script (`S#`).

The legacy per-dataset outputs (`AFDB_v4`, `vDulab`, `vDulab2`, `AF3`,
`DY47073`) are retained as the v1 regression baseline and are disabled by
default.

Segment coordinates are **1-based inclusive**; `length = end - start + 1`.

## Layout

```
config/config.yaml            Path roots, pinned resource versions, parameters
workflow/Snakefile            Config root, stage includes, default target
workflow/rules/<nn>_<stage>.smk    Rules grouped by pipeline stage
workflow/datasets/<nn>_<ds>.smk    Dataset declarations (input paths, version)
workflow/envs/                Driver environment and the shared rule environment
src/                          Functional core + typer CLI
scripts/afdb_domains.py       Thin executable shell
resources/                    PomBase-only curated inputs tracked in git
data/                         Working directory (git-ignored)
results/                      Pipeline outputs (git-ignored)
tests/                        Unit and golden tests
```

## Quickstart

```bash
# Driver environment (Snakemake only)
mamba env create -f workflow/envs/environment.yaml
mamba activate afdb-domains-driver

# Validate the DAG without touching data
snakemake -n --cores 1

# Run the full pipeline (rule environments are created on first use)
snakemake --cores 16 --use-conda
```

Manual CLI use is also supported:

```bash
python scripts/afdb_domains.py --help
python scripts/afdb_domains.py pae-to-domains --pae-dir <dir> --output-dir <dir>
```

## Configuration

`config/config.yaml` declares the machine-specific path roots:

- `cache_root` — PomBase releases and AlphaFold structures.
- `reference_root` — `/hugedata/YushengYang/AFdb`, holding the raw `TED/` and
  `ECOD/` dumps.
- `legacy_root` — the v1 repository, used as the read-only numerical
  verification baseline.
- `raw_root` — the v1 `data/0_raw` directory holding custom predictions.

All other paths are written with a `cache:`, `reference:`, `legacy:`, `raw:`,
or `project:` prefix and resolved at load time.

Resource versions are pinned per role because the legacy results mixed
releases:

| Role | Release | Used for |
|---|---|---|
| PomBase coding | `2026-06-01` | Gene table, PeptideStats, peptide FASTA, GFF3 |
| PomBase features | `2026-06-01` | TED/ECOD protein ID mapping |
| AlphaFold DB | v6 (hugedata) | Canonical structures and PAE |
| AFdb v4 | `2023-07-08` | TED/ECOD coordinate validation only |
| TED | Zenodo 13908086 | Domain assignments |
| ECOD | `dpam_v5` | Domain assignments |

## Pipeline

```
AFdb v6 + PomBase 2026-06-01 ──compare-afdb-pombase──> comparison tables
                                      │
                                      v
                    build-registry ──> gene_structure_registry.tsv
                                      │
        ┌─────────────────────────────┴─────────────────────────────┐
        v                                                           v
registry-pae-to-domains ──> *.domain ──registry-domains-to-table──> PAE table
registry-extract-plddt ──────────────────────────────────────────> pLDDT table
        └────────────── registry-merge-domains ────────────────────┘
                                      │
                                      v
                    registry-unify-segments ──> results/09_segments.tsv
```

The legacy per-dataset path is kept for regression testing:

```
PAE JSON ──pae-to-domains──> *.domain ──domains-to-table──> <ds>_PAE_domain.tsv
AF structures ──extract-plddt──> <ds>_pLDDT_scores.tsv
        └──────── merge-domains ────────┘
                        │
                        v
              <ds>_reformated.tsv ──assign-type──> <ds>_annotated.tsv
                                                        │
TED dump ──ted-pombe──> pombe_ted_segments.tsv ─────────┤
ECOD dump ──ecod-pombe──> pombe_ecod_segments.tsv ──────┴──unify-segments──> segments.tsv
```

The DNA branch (`dna-bed` → `dna-map` → `dna-organize` → `dna-genome-bed`) and
the 3D-AF-Surfer branch (`confident-domains`) are maintained but are not part
of the acceptance criteria.

## Provenance

The legacy repository is tagged `v1-legacy`. Anything not listed under
Deliverables above — the DIT_HAP/CKHS comparison, coding-gene summaries,
custom-prediction FASTA sets, and the earlier PDF figures — is archived there
and intentionally not reproduced by this pipeline. See
`docs/output_manifest.md` for the full inventory.
