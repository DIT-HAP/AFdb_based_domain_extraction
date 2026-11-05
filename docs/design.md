# Design decisions (2026-11-03)

## Version pinning

The v1 results were not produced with a single set of versions. v2 pins each
role instead of pinning one global release:

| Role | Version | Why |
|---|---|---|
| PomBase coding | `2026-06-01` | Gene table, PeptideStats, peptide FASTA, GFF3 |
| PomBase features | `2026-06-01` | Protein ID mapping for the TED/ECOD subsets |
| AlphaFold DB | v6 (hugedata) | Canonical structures and PAE |
| AFdb v4 | `2023-07-08` | Only for TED/ECOD coordinate validation |
| TED | Zenodo 13908086 | Domain assignments |
| ECOD | `dpam_v5` | Domain assignments |

The legacy `AFDB_v4`, `vDulab`, `vDulab2`, `AF3`, and `DY47073` datasets are
declared but disabled; they remain as the v1 regression baseline. The
canonical output comes from the gene structure registry.

## Gene structure registry

The registry selects exactly one structure and one PAE matrix per PomBase
protein-coding gene:

| Rule | Genes |
|---|---:|
| Identical PomBase/AFdb sequence | 5,008 (AFdb v6 primary, including 9 sequence-matched entries) |
| AFdb v6 entry whose sequence equals PomBase | 62 (alternate accession) |
| vDulab prediction whose input equals PomBase | 46 (38 non-identical + 8 long proteins) |
| AF3 prediction (long proteins without vDulab) | 3 (dhc1, tra2, mdn1) |
| Excluded | 7 (SPCC569.03, sib1, B5BP45-49) |
| AFdb-only (annotated, not analysed) | 16 |

`in_analysis` is true for 5,119 genes. The registry also records
`structure_sequence`, `modeled_length`, `coordinate_basis`, and
`ted_ecod_coordinate_ok`. The TED/ECOD flag is exact sequence equality between
the AFdb v4 model and the selected structure; genes with a v4/v6 difference
are excluded from TED/ECOD rather than risking coordinate drift.

Python dependencies are pinned in `workflow/envs/core.yaml`. The critical pin
is `python-igraph=0.11.4`: version 1.0 produces a different Leiden partition
(AF3 yielded 181 domains instead of 171). `pandas=2.2.3`, `numpy=1.26.4`, and
`biopython=1.85` are the versions that reproduced the legacy tables exactly.

## Workflow layout

Rules are grouped by pipeline stage in `workflow/rules/`: `01_registry.smk`
(canonical chain), `02_pae.smk`, `03_plddt.smk`, `04_domains.smk`,
`05_external.smk`, `06_unify.smk`, and `07_dna.smk`. The root `Snakefile` only
holds configuration, includes, and `rule all`. Dataset-specific input paths
and version labels live at the top of `02_pae.smk`, the entry point of the
legacy per-dataset chain (D1–D5 in `docs/results_index.md`). The unused
`AFDB_20251107` declaration was dropped; the canonical registry reads the v6
structures through `config/config.yaml`.

## Resource policy

- `resources/` holds small PomBase-only inputs that the pipeline reads directly
  and that are committed: PomBase tables, the TED/ECOD PomBase subsets, the
  DY47073 protein subset, and the provenance manifest.
- Raw upstream dumps live under the reference root
  (`/hugedata/YushengYang/AFdb/TED` and `ECOD`): the 128G TED domain summary,
  the 27G TED per-chain segmentation, and the 31M ECOD dump. AlphaFold
  structures and PAE stay in the cache. `config/config.yaml` points at them;
  nothing large is committed.
- `ted-pombe` and `ecod-pombe` stream the raw dumps and write the PomBase
  subsets into `resources/`, so the extraction is reproducible from the
  workflow alone.
- `fetch-resources` records size and SHA-256 (files under 50MB) for every
  configured resource.

## Unified segment schema

One long-format table per dataset, plus a combined table:

```
dataset, uniprot_id, systematic_id, gene_name, protein_length,
source, source_version, segment_id, segment_part,
start, end, length, score, score_type, fraction_gt_70,
segment_class, family, location
```

- `uniprot_id` is the only join key across the three sources.
- Coordinates are 1-based inclusive; `length = end - start + 1`.
- Multi-interval domains are exploded; `segment_part` numbers the parts.
- Missing values are empty strings in the TSV.
- The three sources are concatenated, not inner-joined, because only 4,216 of
  the 5,128 AFDB proteins have both TED and ECOD assignments.
- TED and ECOD are filtered to the proteins predicted by the dataset. A
  vDulab2 table (71 proteins) therefore has 700 rows, not the 22k rows of a
  full-proteome table.

## Scope

- The DIT_HAP and CKHS analyses are out of scope. Their 16 PDF figures and the
  notebook that produced them are archived in `v1-legacy`.
- The exploratory notebooks are archived. Their deterministic logic moved into
  CLI subcommands and Snakemake rules.
- The DNA branch and the 3D-AF-Surfer confident-domain branch are maintained
  and runnable, but they are not part of the acceptance criteria.

## Verification

| Dataset | Result |
|---|---|
| AF3 | Exact match on all four PAE tables |
| vDulab2 | Exact on intervals, lengths, fractions, types; 2 rows differ by ±0.01 mean_pLDDT |
| vDulab | Exact on the same fields; 3 rows differ by ±0.01 mean_pLDDT (after pinning the 2024-03-05 input release) |
| AFDB v4 | PAE and pLDDT tables exact; 17,851 unique domains after removing 9 duplicate legacy rows; 68 rows differ by ±0.01 mean_pLDDT; 42 rows differ only by the intended UniProt fallback in `Systematic ID` |
| confident-domains | 20,554 unique rows match the legacy table; the legacy file duplicated 3 proteins and drifted by ±0.01 on 45 rows |

The `Systematic ID` fallback is present in the current v1 `merge_PAE_domains.py`
but missing from the 2023 output file. v2 follows the current script: unmapped
proteins keep their UniProt accession instead of an empty systematic id.
