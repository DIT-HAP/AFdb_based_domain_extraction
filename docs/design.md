# Design decisions (2026-11-03)

## Version pinning

The v1 results were not produced with a single set of versions. v2 pins each
role instead of pinning one global release:

| Role | Version | Why |
|---|---|---|
| PomBase coding | `2023-07-01` | Gene table, PeptideStats, GFF3 used by the PAE branch |
| PomBase features | `2025-05-01` | Protein ID mapping for the TED/ECOD subsets |
| AlphaFold DB | `2023-07-08` | AFDB v4 structures and PAE |
| TED | Zenodo 13908086 | Domain assignments |
| ECOD | `dpam_v5` | Domain assignments |

Python dependencies are pinned in `workflow/envs/core.yaml`. The critical pin
is `python-igraph=0.11.4`: version 1.0 produces a different Leiden partition
(AF3 yielded 181 domains instead of 171). `pandas=2.2.3`, `numpy=1.26.4`, and
`biopython=1.85` are the versions that reproduced the legacy tables exactly.

## Workflow layout

Rules are grouped by pipeline stage in `workflow/rules/`: `pae.smk`,
`plddt.smk`, `domains.smk`, `external.smk`, `unify.smk`, and `dna.smk`. The
root `Snakefile` only holds configuration, dataset declarations, includes, and
`rule all`. Dataset-specific input paths and version labels live in
`workflow/datasets/<dataset>.smk`.

## Resource policy

- `resources/` holds small inputs that the pipeline reads directly and that are
  committed: PomBase tables, the curated TED PomBase table, the TED/ECOD PomBase
  subsets, the DY47073 protein subset, and the provenance manifest.
- Large upstream dumps stay in the cache: TED 146G, ECOD 30M, AlphaFold
  structures and PAE. `config/config.yaml` points at them; nothing large is
  copied or committed.
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
