# Results index

Numbered inventory of the active v2 outputs under `results/`, in pipeline
order. Row counts exclude the header line and reflect the 2026-11-05 run.

Four numbering schemes are used:

- `#n` — output file or directory; the number is also the on-disk filename
  prefix (`01_afdb/`, `04_gene_structure_registry.tsv`, ...).
- `Rn` — Snakemake rule (Appendix A); the rule files are numbered
  `workflow/rules/01_afdb.smk`–`11_dna.smk`.
- `Sn` — script or Python module (Appendix B).
- `Dn` — dataset declaration (Appendix C), defined in the `DATASETS` block at
  the top of `workflow/rules/06_pae.smk`.

Canonical deliverables: **#4** (`04_gene_structure_registry.tsv`) and **#11**
(`11_segments.tsv`). **#1** (`03_afdb_entries.tsv`) is the reusable AFdb
inventory. Everything else is an intermediate or a diagnostic.

## Input branch

### 1. AFdb structure inventory

Producers: **R1** `afdb_cif_entries`, **R2** `afdb_pae_index`,
**R3** `afdb_entry_table` → **S1** `scripts/afdb_domains.py` → **S3**
`src/cli.py` → **S12** `src/afdb.py` → `results/01_afdb/`

| # | File | Rows | Contents |
|---|---|---:|---|
| 1.1 | `01_afdb_cif_entries.tsv` | 5,196 | Accession, description, modelled length/sequence, CIF filename |
| 1.2 | `02_afdb_pae_index.tsv` | 5,196 | Accession and PAE filename |
| 1.3 | `03_afdb_entries.tsv` | 5,196 | 1.1 left-joined with 1.2; 5,196 entries have PAE |

PAE files without a CIF entry are logged by `afdb_entry_table` and excluded
from 1.3.

### 2. UniProt refresh candidates (auxiliary)

Producers: **R4** `refresh_uniprot_triangulation`, **R5**
`refresh_uniprot_resolution` → **S1** → **S3** → **S13** `src/uniprot.py` →
`results/02_uniprot/`

These rules are not in `rule all`. They read `03_comparison/` outputs and write
candidates for manual review; the canonical pipeline reads the pinned snapshots
in `resources/uniprot/`.

| # | File | Contents |
|---|---|---|
| 2.1 | `uniprot_triangulation_refresh.tsv` | 4 columns for the non-identical pairs |
| 2.2 | `uniprot_resolution_refresh.tsv` | 7 columns for the AFdb-only accessions |

## Registry branch

### 3. AFdb vs PomBase comparison

Producer: **R6** `compare_afdb_pombase` → **S2** `scripts/compare_afdb_pombase.py`
→ `results/03_comparison/`

The `01_`–`10_` file prefixes match rows 3.1–3.10 below.

| # | File | Rows | Contents |
|---|---|---:|---|
| 3.1 | `01_merged_table.tsv` | 5,142 | Per-gene AFdb/PomBase/UniProt merge; input to the registry |
| 3.2 | `02_afdb_entry_roles.tsv` | 5,196 | AFdb v6 entry roles against protein entities |
| 3.3 | `03_shared_comparison.tsv` | 5,100 | Genes present in both sources, field by field |
| 3.4 | `04_shared_by_sequence.tsv` | 9 | Entries linked by exact sequence only |
| 3.5 | `05_nonidentical_pairs.tsv` | 101 | Pairs whose sequences differ |
| 3.6 | `06_afdb_only.tsv` | 100 | Entries without a PomBase gene |
| 3.7 | `07_pombase_only.tsv` | 26 | Genes without an AFdb structure |
| 3.8 | `08_summary_counts.tsv` | 45 | Headline counts |
| 3.9 | `09_report.md` | — | Readable report with numbered sections |
| 3.10 | `10_afdb_vs_pombase.xlsx` | — | Excel workbook of the same tables |

### 4. Gene structure registry

Producer: **R7** `build_registry` → **S1** → **S3** → **S5** `src/registry.py`
→ `results/04_gene_structure_registry.tsv`

5,142 rows. One row per gene, recording the selected structure, PAE matrix,
modelled sequence, coordinate basis, and `in_analysis` flag.

### 5. Registry summary

Producer: **R7** (same rule as #4) → `results/05_registry_summary.tsv`

9 rows. Gene counts grouped by `selection_reason` and `structure_source`.

## Curated branch

Every step below reads the registry (#4) and feeds the next step.

### 6. Raw PAE domains

Producer: **R8** `registry_pae_to_domains` → **S1** → **S3** → **S6** `src/pae.py`
→ `results/06_domains/`

5,120 `.domain` files, one per analysed gene.

### 7. PAE domain table

Producer: **R9** `registry_domains_to_table` → **S1** → **S3** → **S6**
→ `results/07_PAE_domain.tsv`

30,757 rows. The `.domain` files (#6) collapsed into one table.

### 8. pLDDT table

Producer: **R10** `registry_extract_plddt` → **S1** → **S3** → **S7**
`src/plddt.py` → `results/08_pLDDT_scores.tsv`

5,119 rows. Per-gene, per-residue confidence scores.

### 9. Merged PAE + pLDDT table

Producer: **R11** `registry_merge_domains` → **S1** → **S3** → **S8**
`src/segments.py` → `results/09_reformated_PAE_domain.tsv`

17,843 rows. Domain intervals joined to their pLDDT statistics.

### 10. Annotated domain table

Producer: **R12** `registry_assign_type` → **S1** → **S3** → **S8**
→ `results/10_annotated_PAE_domain.tsv`

17,843 rows. Adds the segment type (`disordered`, `structured`, and so on).

### 11. Unified segment table

Producer: **R13** `registry_unify_segments` → **S1** → **S3** → **S9**
`src/unify.py` → `results/11_segments.tsv`

42,178 rows: PAE 20,711 + TED 11,983 + ECOD 9,484. This is the final delivery
table; coordinates are 1-based inclusive.

## Archived

### 12. Stale 2026-11-03 dataset outputs

Location: `tmp/afdb_vs_pombase/archive/stale_results_20261103/`

Holds `DY47073`, `vDulab2`, `vDulab`, `AF3`, `AFDB_v4`, and `combined`.
These datasets are declared with `enabled: False` and are not part of the
canonical run. They were moved out of `results/` on 2026-11-05 and can be
regenerated by enabling the corresponding declaration.

The pre-split canonical outputs (`results/01_comparison/`,
`02_gene_structure_registry.tsv`, `03_registry_summary.tsv`,
`09_segments.tsv`) remain on disk as the byte-equality baseline until the
full chain is rerun.

## Appendix A — Rule index

`R1`–`R13` form the active canonical chain. `R4`/`R5` are auxiliary refresh
rules outside `rule all`. `R14`–`R23` are the legacy per-dataset path
(`enabled: False`). `R24`–`R27` are the DNA branch, which is maintained but
outside the acceptance criteria.

| # | Rule | Snakefile | Entry point | Script(s) | Produces |
|---|---|---|---|---|---|
| R1 | `afdb_cif_entries` | `rules/01_afdb.smk` | `afdb-cif-entries` | S1, S3, S12 | #1.1 |
| R2 | `afdb_pae_index` | `rules/01_afdb.smk` | `afdb-pae-index` | S1, S3, S12 | #1.2 |
| R3 | `afdb_entry_table` | `rules/01_afdb.smk` | `afdb-entry-table` | S1, S3, S12 | #1.3 |
| R4 | `refresh_uniprot_triangulation` | `rules/02_uniprot.smk` | `refresh-uniprot-triangulation` | S1, S3, S13 | #2.1 |
| R5 | `refresh_uniprot_resolution` | `rules/02_uniprot.smk` | `refresh-uniprot-resolution` | S1, S3, S13 | #2.2 |
| R6 | `compare_afdb_pombase` | `rules/03_compare.smk` | `compare_afdb_pombase.py` | S2 | #3 |
| R7 | `build_registry` | `rules/04_registry.smk` | `build-registry` | S1, S3, S5 | #4, #5 |
| R8 | `registry_pae_to_domains` | `rules/05_curated.smk` | `registry-pae-to-domains` | S1, S3, S6 | #6 |
| R9 | `registry_domains_to_table` | `rules/05_curated.smk` | `registry-domains-to-table` | S1, S3, S6 | #7 |
| R10 | `registry_extract_plddt` | `rules/05_curated.smk` | `registry-extract-plddt` | S1, S3, S7 | #8 |
| R11 | `registry_merge_domains` | `rules/05_curated.smk` | `registry-merge-domains` | S1, S3, S8 | #9 |
| R12 | `registry_assign_type` | `rules/05_curated.smk` | `assign-type` | S1, S3, S8 | #10 |
| R13 | `registry_unify_segments` | `rules/05_curated.smk` | `registry-unify-segments` | S1, S3, S9 | #11 |
| R14 | `pae_to_domains` | `rules/06_pae.smk` | `pae-to-domains` | S1, S3, S6 | `<ds>/domains/` |
| R15 | `domains_to_table` | `rules/06_pae.smk` | `domains-to-table` | S1, S3, S6 | `<ds>/PAE_domain.tsv` |
| R16 | `extract_plddt` | `rules/07_plddt.smk` | `extract-plddt` | S1, S3, S7 | `<ds>/pLDDT_scores.tsv` |
| R17 | `merge_domains` | `rules/08_domains.smk` | `merge-domains` | S1, S3, S8 | `<ds>/reformated_PAE_domain.tsv` |
| R18 | `assign_type` | `rules/08_domains.smk` | `assign-type` | S1, S3, S8 | `<ds>/annotated_PAE_domain.tsv` |
| R19 | `confident_domains` | `rules/08_domains.smk` | `confident-domains` | S1, S3, S8 | `<ds>/confident_domains.tsv` |
| R20 | `ted_pombe` | `rules/09_external.smk` | `ted-pombe` | S1, S3, S10 | `resources/ted/pombe_ted_segments.tsv` |
| R21 | `ecod_pombe` | `rules/09_external.smk` | `ecod-pombe` | S1, S3, S10 | `resources/ecod/pombe_ecod_segments.tsv` |
| R22 | `unify_segments` | `rules/10_unify.smk` | `unify-segments` | S1, S3, S9 | `<ds>/segments.tsv` |
| R23 | `combine_segments` | `rules/10_unify.smk` | `combine-segments` | S1, S3, S9 | `combined/segments.tsv` |
| R24 | `dna_bed` | `rules/11_dna.smk` | `dna-bed` | S1, S3, S11 | `data/nucleotide_bed.tsv` |
| R25 | `dna_map` | `rules/11_dna.smk` | `dna-map` | S1, S3, S11 | `<ds>/DNA_domain.tsv` |
| R26 | `dna_organize` | `rules/11_dna.smk` | `dna-organize` | S1, S3, S11 | `<ds>/merged_DNA_domain.tsv` |
| R27 | `dna_genome_bed` | `rules/11_dna.smk` | `dna-genome-bed` | S1, S3, S11 | `<ds>/domains.bed` |

## Appendix B — Script index

| # | Script | Role | Used by |
|---|---|---|---|
| S1 | `scripts/afdb_domains.py` | Executable shell; adds the repo root to `sys.path` and calls `src.cli.main` | R1–R5, R7–R27 |
| S2 | `scripts/compare_afdb_pombase.py` | Standalone AFdb-vs-PomBase comparison; writes tables, report, and workbook | R6 |
| S3 | `src/cli.py` | Typer application; one subcommand per stage, exit-code handling | R1–R5, R7–R27 |
| S4 | `src/files.py` | Config loading, logging, TSV read/write, required-column checks | All CLI commands |
| S5 | `src/registry.py` | Structure selection, sequence matching, registry assembly | R7 |
| S6 | `src/pae.py` | PAE parsing/clustering, `.domain` files, domain table | R8, R9, R14, R15 |
| S7 | `src/plddt.py` | pLDDT extraction, interval utilities, domain statistics | R10, R16 |
| S8 | `src/segments.py` | Domain merge, type annotation, confident-domain caller | R11, R12, R17–R19 |
| S9 | `src/unify.py` | Segment loading, PomBase annotation, per-dataset and combined tables | R13, R22, R23 |
| S10 | `src/external.py` | PomBase subsets of the TED and ECOD dumps | R20, R21 |
| S11 | `src/dna.py` | GFF parsing, nucleotide projection, merged BED | R24–R27 |
| S12 | `src/afdb.py` | CIF/PAE scanning and the AFdb entry inventory | R1–R3 |
| S13 | `src/uniprot.py` | UniProt REST fetch and refresh assembly | R4, R5 |

## Appendix C — Dataset declaration index

Declarations live in the `DATASETS` block at the top of
`workflow/rules/06_pae.smk`. All five are currently `enabled: False`; the
canonical run does not read them. The planned `AFDB_20251107` (AFdb v6) entry
was dropped: nothing referenced it, and the registry reads v6 through
`config/config.yaml`.

| # | Dataset | Source version |
|---|---|---|
| D1 | `AFDB_v4` | `AFDB_v4_20230708` |
| D2 | `vDulab` | `vDulab_20240305` |
| D3 | `vDulab2` | `vDulab2_20240904` |
| D4 | `AF3` | `AF3_20241101` |
| D5 | `DY47073` | `vDulab2_20240904` |
