# AFdb-based domain extraction

Call protein domains from AlphaFold PAE matrices, merge them with per-residue
pLDDT confidence, and emit one unified segment table per dataset alongside the
PomBase subsets of TED and ECOD.

The v2 refactor keeps the original numerical behaviour where it was defined and
replaces the ad-hoc shell orchestration with a Snakemake workflow.

## Deliverables

- `results/01_afdb/03_afdb_entries.tsv` — one row per AFdb v6 entry with the
  modelled sequence, description, CIF filename, and PAE availability.
- `results/04_registry/01_gene_structure_registry.tsv` — one row per PomBase coding
  gene (plus the AFdb-only entries), recording the selected structure, PAE
  matrix, modelled sequence, coordinate basis, and analysis flag.
- `results/04_registry/03_registry_stats.tsv` — downstream counts derived from
  the registry (LocalAF2_PomBaseDiffAFdb, AF3_LongProt, and TED/ECOD coverage).
- `results/10_segments.tsv` — the canonical unified segment table with
  three sources: `PAE`, `TED`, `ECOD`.
- `results/03_comparison/` — the detailed AFdb-vs-PomBase tables, report,
  and Excel workbook.

`docs/results_index.md` numbers every output file, its row count, the producing
rule (`R#`), and the script (`S#`).

The legacy per-dataset branch was removed on 2026-11-19; its archived outputs
remain under `tmp/afdb_vs_pombase/archive/` and the v1 code in git history.

Segment coordinates are **1-based inclusive**; `length = end - start + 1`.

## Layout

```
config/config.yaml            Path roots, pinned resource versions, parameters
workflow/Snakefile            Config root, shared result paths, stage includes
workflow/rules/01_afdb.smk    AFdb structure inventory
workflow/rules/02_uniprot.smk UniProt snapshot refresh (manual, not in rule all)
workflow/rules/03_compare.smk AFdb vs PomBase comparison
workflow/rules/04_registry.smk Gene structure registry
workflow/rules/05_curated.smk PAE clustering through unified segments
workflow/rules/06_external.smk TED/ECOD PomBase subsets
workflow/rules/07_dna.smk      DNA projection branch
workflow/envs/                Driver environment and the shared rule environment
src/                          Functional core + typer CLI
scripts/afdb_domains.py       Thin executable shell
scripts/compare_afdb_pombase.py Standalone comparison script
resources/                    PomBase-only curated inputs tracked in git
data/                         Working directory (git-ignored)
results/                      Pipeline outputs (git-ignored)
tests/                        Unit tests
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
python scripts/afdb_domains.py registry-pae-to-domains --registry <tsv> --output-dir <dir>
```

## Configuration

`config/config.yaml` declares the machine-specific path roots:

- `cache_root` — local scratch root, unused by the pipeline.
- `reference_root` — `/hugedata/YushengYang/AFdb`, holding the AFdb structures
  and PAE, the selected custom predictions under `predictions/`, and the raw
  `TED/` and `ECOD/` dumps.
- `legacy_root` — the v1 repository, used as the read-only numerical
  verification baseline.

All other paths are written with a `cache:`, `reference:`, `legacy:`, or
`project:` prefix and resolved at load time.

Resource versions are pinned per role because the legacy results mixed
releases:

| Role | Release | Used for |
|---|---|---|
| PomBase coding | `2026-06-01` | Gene table, PeptideStats, peptide FASTA, GFF3 |
| PomBase features | `2026-06-01` | TED/ECOD protein ID mapping |
| AlphaFold DB | v6 (hugedata) | Canonical structures and PAE |
| AFdb v4 | `2023-07-08` | TED/ECOD coordinate validation only |
| Local AF2 predictions | `2024-03-05` | 45 models for PomBase-vs-AFdb discordant genes |
| AlphaFold 3 | `2024-11-01` | 3 long proteins without a usable local AF2 model |
| TED | Zenodo 13908086 | Domain assignments |
| ECOD | `dpam_v5` | Domain assignments |

Selected custom predictions live under `reference_root/predictions/<version>/`
with `structure/`, `pae/`, and `input_fasta/` subdirectories; the top-level
`manifest.tsv` records source files and per-gene selection status.

## Pipeline

```
AFdb v6 CIF + PAE ──afdb-cif-entries / afdb-pae-index──> AFdb entry inventory
                                                              │
PomBase 2026-06-01 + UniProt snapshots ───────────────────────┤
                                                              v
                                              compare-afdb-pombase ──> comparison tables
                                                              │
                                                              v
                            build-registry ──> results/04_registry/01_gene_structure_registry.tsv
                                                              │
        ┌─────────────────────────────────────────────────────┴─────┐
        v                                                           v
registry-pae-to-domains ──> *.domain ──registry-domains-to-table──> PAE table
registry-extract-plddt ──────────────────────────────────────────> pLDDT table
        └────────────── registry-merge-domains ────────────────────┘
                                      │
                                      v
                    registry-unify-segments ──> results/10_segments.tsv
```

The two UniProt snapshots are pinned resources. `02_uniprot.smk` provides
manual refresh rules that write candidates to `results/02_uniprot/`; they are
not part of `rule all` (see `resources/README.md`).

The DNA branch (`dna-bed` → `dna-map` → `dna-organize` → `dna-genome-bed`)
reads the canonical annotated table and writes `results/12_dna/`. The
3D-AF-Surfer branch (`registry-confident-domains`) reads the canonical pLDDT
table plus the registry and writes `results/11_confident_domains.tsv`. Both are
maintained but are not part of the acceptance criteria.

## Provenance

The legacy repository is tagged `v1-legacy`. Anything not listed under
Deliverables above — the DIT_HAP/CKHS comparison, coding-gene summaries,
custom-prediction FASTA sets, and the earlier PDF figures — is archived there
and intentionally not reproduced by this pipeline. See
`docs/output_manifest.md` for the full inventory.
