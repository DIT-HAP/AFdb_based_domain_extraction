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
protein-coding gene. The comparison stage pairs genes and AFdb entries by
protein sequence first and falls back to the UniProt accession only when no
exact sequence match exists:

| Rule | Genes |
|---|---:|
| Sequence-identical AFdb v6 entry | 5,071 (4,999 own accession + 63 alternate + 9 without an AFdb entry at their accession) |
| vDulab prediction whose input equals PomBase | 45 (37 non-identical + 8 long proteins) |
| AF3 prediction (long proteins without vDulab) | 3 (dhc1, tra2, mdn1) |
| Excluded | 7 genes (SPCC569.03, sib1, B5BP45-49) + 16 AFdb-only rows |
| Unselected AFdb entries | 92 (63 superseded by a sequence match, 13 duplicates, 16 no PomBase match) |

`in_analysis` is true for 5,119 genes. The registry also records
`structure_sequence`, `modeled_length`, `coordinate_basis`, and
`ted_ecod_coordinate_ok`. The TED/ECOD flag is exact sequence equality between
the AFdb v4 model and the selected structure; genes with a v4/v6 difference
are excluded from TED/ECOD rather than risking coordinate drift.

Sequence-first pairing changes the source of 63 genes whose AFdb v6 primary
entry no longer matches the PomBase peptide: the registry uses the alternate
entry whose sequence is exact. One of these is hht3 (SPBC1105.11c), whose
PomBase 2026 peptide equals the H3.1/H3.2 sequence (P09988) rather than the
current UniProt H3.3 entry (P10651); the registry therefore selects P09988 and
marks the coordinate basis as PomBase.

Python dependencies are pinned in `workflow/envs/core.yaml`. The critical pin
is `python-igraph=0.11.4`: version 1.0 produces a different Leiden partition
(AF3 yielded 181 domains instead of 171). `pandas=2.2.3`, `numpy=1.26.4`, and
`biopython=1.85` are the versions that reproduced the legacy tables exactly.

## Workflow layout

Rules are grouped by pipeline stage, one file per stage, in `workflow/rules/`:

| File | Stage |
|---|---|
| `01_afdb.smk` | AFdb CIF/PAE inventory (`afdb_cif_entries`, `afdb_pae_index`, `afdb_entry_table`) |
| `02_uniprot.smk` | Manual UniProt snapshot refresh (not in `rule all`) |
| `03_compare.smk` | AFdb vs PomBase comparison |
| `04_registry.smk` | Gene structure registry |
| `05_curated.smk` | PAE clustering through unified segments |
| `06_pae.smk` | Legacy dataset declarations and PAE rules (D1–D5) |
| `07_plddt.smk`–`11_dna.smk` | Legacy pLDDT, domain, external, unification, DNA branches |

The root `Snakefile` holds the configuration, the shared numbered result path
constants, the stage includes, and `rule all`. Dataset-specific input paths
and version labels live at the top of `06_pae.smk`, the entry point of the
legacy per-dataset chain. The unused `AFDB_20251107` declaration was dropped;
the canonical registry reads the v6 structures through `config/config.yaml`.

The AFdb inventory is produced once and cached, so the comparison and registry
no longer parse 5,196 gzip CIF headers themselves. `03_afdb_entries.tsv`
carries `afdb_accession`, `afdb_description`, `afdb_length`, `afdb_sequence`,
`cif_file`, `pae_available`, and `pae_file`. The join starts from the CIF
inventory; PAE files without a CIF entry are reported and ignored.

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
- The two UniProt tables stay committed snapshots. `02_uniprot.smk` adds
  `refresh_uniprot_triangulation` and `refresh_uniprot_resolution`; both are
  manual rules outside `rule all`, write minimal-column candidates to
  `results/02_uniprot/`, and require a prior comparison run for their query
  list. Promotion (copy to `resources/uniprot/`, update `config/config.yaml`)
  is documented in `resources/README.md`.

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

## Stage split verification (2026-11-05)

`01_registry.smk` was split into `01_afdb.smk`–`05_curated.smk`, and the
comparison script stopped parsing CIF/PAE directories. The front end was rerun
from `results/01_afdb/03_afdb_entries.tsv`: all seven comparison TSVs, the
entry-role table, the registry, and the summary are byte-identical to the
pre-split outputs; the Excel workbook has the same 15 sheets with identical
values. The Markdown report differs only in the corrected registry path
reference (`04_gene_structure_registry.tsv`). The PAE clustering/downstream
chain was not rerun for this check.

The two refresh rules were also run end to end. Their candidates reproduce the
pinned snapshots exactly on every column the pipeline reads: 100/100 resolution
rows and 101/101 triangulation rows.

## Sequence-first matching (2026-11-05)

The comparison now pairs genes to AFdb entries by exact protein sequence before
falling back to the UniProt accession. The comparison and registry were rerun:

- 5,109 paired genes: 5,071 sequence matches + 38 UniProt ID fallbacks.
- 72 sequence matches use an entry other than the gene's UniProt accession
  (63 superseded old entries + 9 genes without an AFdb entry at their accession).
- 92 unselected AFdb entries: 63 superseded, 13 duplicates, 16 no PomBase match.
- 17 PomBase-only genes; merged table still 5,142 rows.
- Registry `in_analysis` remains 5,119 (5,071 AFdb + 45 vDulab + 3 AF3); hht3
  moves from vDulab to the sequence-identical AFdb entry P09988.
- The Excel workbook gains an `afdb_only_superseded` sheet.

The earlier byte-identical stage-split check was performed before this matching
change and no longer applies verbatim.
