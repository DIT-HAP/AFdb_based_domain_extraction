# Results index

Numbered inventory of the active v2 outputs under `results/`, in pipeline
order. Row counts exclude the header line and reflect the 2026-11-05 run.

Four numbering schemes are used:

- `#n` — output file or directory; the number is also the on-disk filename
  prefix (`01_afdb/`, `04_registry/`, ...).
- `Rn` — Snakemake rule (Appendix A); rule files are numbered
  `workflow/rules/01_afdb.smk`–`11_dna.smk` (06–08 and 10 were removed).
- `Sn` — script or Python module (Appendix B).
- `Dn` — dataset declaration (Appendix C), defined in the `datasets` block of
  `config/config.yaml` and loaded by `src/files.py`.

> 2026-11-19: the legacy per-dataset branch (`06_pae.smk`, `07_plddt.smk`,
> `08_domains.smk`, `10_unify.smk`) was removed. The archived outputs below and
> the v1 code in git history remain the only way to reproduce those tables.

Canonical deliverables: **#4** (`04_registry/01_gene_structure_registry.tsv`)
and **#10** (`10_segments.tsv`). **#1** (`03_afdb_entries.tsv`) is the
reusable AFdb inventory. Everything else is an intermediate or a diagnostic.

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

The `01_`–`12_` file prefixes match rows 3.1–3.12 below.

| # | File | Rows | Contents |
|---|---|---:|---|
| 3.1 | `01_merged_table.tsv` | 5,142 | Per-gene sequence-first AFdb/PomBase merge with both sequences and `relationship`; input to the registry |
| 3.2 | `02_afdb_entry_roles.tsv` | 5,196 | AFdb v6 entry roles against selected genes |
| 3.3 | `03_shared_comparison.tsv` | 5,109 | Paired genes, field by field (5,071 sequence + 38 ID fallback) |
| 3.4 | `04_shared_by_sequence.tsv` | 72 | Sequence matches whose entry is not the gene's UniProt accession |
| 3.5 | `05_nonidentical_pairs.tsv` | 38 | ID-fallback pairs whose sequences differ |
| 3.6 | `06_afdb_only.tsv` | 92 | AFdb entries selected by no gene |
| 3.7 | `07_pombase_only.tsv` | 17 | Genes without an AFdb structure |
| 3.8 | `08_summary_counts.tsv` | 42 | Headline counts |
| 3.9 | `09_report.md` | — | Readable report with numbered sections |
| 3.10 | `10_afdb_vs_pombase.xlsx` | — | Excel workbook of the same tables |
| 3.11 | `11_pombase_coding_genes.tsv` | 5,126 | Protein-coding gene table used for pairing |
| 3.12 | `12_pombase_peptides.tsv` | 5,126 | `peptide.fa` converted to gene ID and raw sequence (stops retained) |

### 4. Gene structure registry

Producer: **R7** `build_registry` → **S1** → **S3** → **S5** `src/registry.py`
→ `results/04_registry/`

The same rule writes all three files.

| # | File | Rows | Contents |
|---|---|---:|---|
| 4.1 | `01_gene_structure_registry.tsv` | 5,142 | One row per gene: selected structure, PAE matrix, modelled sequence, coordinate basis, and `in_analysis` flag |
| 4.2 | `02_registry_summary.tsv` | 9 | Gene counts grouped by `selection_reason` and `structure_source` |
| 4.3 | `03_registry_stats.tsv` | 8 | Downstream counts derived from the registry (LocalAF2_PomBaseDiffAFdb, AF3_LongProt, TED/ECOD) |

## Curated branch

Every step below reads the registry (#4) and feeds the next step.

### 5. Raw PAE domains

Producer: **R8** `registry_pae_to_domains` → **S1** → **S3** → **S6** `src/pae.py`
→ `results/05_domains/`

5,119 `.domain` files, one per analysed gene.

### 6. PAE domain table

Producer: **R9** `registry_domains_to_table` → **S1** → **S3** → **S6**
→ `results/06_PAE_domain.tsv`

30,757 rows. The `.domain` files (#5) collapsed into one table.

### 7. pLDDT table

Producer: **R10** `registry_extract_plddt` → **S1** → **S3** → **S7**
`src/plddt.py` → `results/07_pLDDT_scores.tsv`

5,119 rows. Per-gene, per-residue confidence scores.

### 8. Merged PAE + pLDDT table

Producer: **R11** `registry_merge_domains` → **S1** → **S3** → **S8**
`src/segments.py` → `results/08_reformated_PAE_domain.tsv`

17,843 rows. Domain intervals joined to their pLDDT statistics.

### 9. Annotated domain table

Producer: **R12** `registry_assign_type` → **S1** → **S3** → **S8**
→ `results/09_annotated_PAE_domain.tsv`

17,843 rows. Adds the segment type (`disordered`, `structured`, and so on).

### 10. Unified segment table

Producer: **R13** `registry_unify_segments` → **S1** → **S3** → **S9**
`src/unify.py` → `results/10_segments.tsv`

42,178 rows: PAE 20,711 + TED 11,983 + ECOD 9,484. This is the final delivery
table; coordinates are 1-based inclusive.

## Archived

### 11. Stale 2026-11-03 dataset outputs

Location: `tmp/afdb_vs_pombase/archive/stale_results_20261103/`

Holds the legacy directories `DY47073`, `vDulab2`, `vDulab`, `AF3`, `AFDB_v4`,
and `combined` (old dataset names). They were moved out of `results/` on
2026-11-05. The rules that produced them were removed on 2026-11-19, so they
can only be regenerated from git history.

### 12. Pre-sequence-first baseline

Location: `results/arc/`

Holds the complete pre-sequence-first numbering: `01_comparison/`,
`02_gene_structure_registry.tsv`, `03_registry_summary.tsv`, and the curated
outputs `04_domains/`, `05_PAE_domain.tsv`, `06_pLDDT_scores.tsv`,
`07_reformated_PAE_domain.tsv`, `08_annotated_PAE_domain.tsv`, and
`09_segments.tsv`. The comparison and registry files were moved out of
`results/` on 2026-11-11; the curated outputs, which derive from that earlier
registry, were moved on 2026-11-16. All are kept for reference only; the
active outputs are #3 and #4. Rerun the curated branch to produce #5–#10
under the current numbering.

## Appendix A — Rule index

`R1`–`R13` form the active canonical chain. `R4`/`R5` are auxiliary refresh
rules outside `rule all`; `R14` is the confident-domain branch, also outside
`rule all`. `R17`–`R20` are the DNA branch, which is maintained but outside the
acceptance criteria.

| # | Rule | Snakefile | Entry point | Script(s) | Produces |
|---|---|---|---|---|---|
| R1 | `afdb_cif_entries` | `rules/01_afdb.smk` | `afdb-cif-entries` | S1, S3, S12 | #1.1 |
| R2 | `afdb_pae_index` | `rules/01_afdb.smk` | `afdb-pae-index` | S1, S3, S12 | #1.2 |
| R3 | `afdb_entry_table` | `rules/01_afdb.smk` | `afdb-entry-table` | S1, S3, S12 | #1.3 |
| R4 | `refresh_uniprot_triangulation` | `rules/02_uniprot.smk` | `refresh-uniprot-triangulation` | S1, S3, S13 | #2.1 |
| R5 | `refresh_uniprot_resolution` | `rules/02_uniprot.smk` | `refresh-uniprot-resolution` | S1, S3, S13 | #2.2 |
| R6 | `compare_afdb_pombase` | `rules/03_compare.smk` | `compare_afdb_pombase.py` | S2 | #3 |
| R7 | `build_registry` | `rules/04_registry.smk` | `build-registry` | S1, S3, S5 | #4.1, #4.2, #4.3 |
| R8 | `registry_pae_to_domains` | `rules/05_curated.smk` | `registry-pae-to-domains` | S1, S3, S6 | #5 |
| R9 | `registry_domains_to_table` | `rules/05_curated.smk` | `registry-domains-to-table` | S1, S3, S6 | #6 |
| R10 | `registry_extract_plddt` | `rules/05_curated.smk` | `registry-extract-plddt` | S1, S3, S7 | #7 |
| R11 | `registry_merge_domains` | `rules/05_curated.smk` | `registry-merge-domains` | S1, S3, S8 | #8 |
| R12 | `registry_assign_type` | `rules/05_curated.smk` | `assign-type` | S1, S3, S8 | #9 |
| R13 | `registry_unify_segments` | `rules/05_curated.smk` | `registry-unify-segments` | S1, S3, S9 | #10 |
| R14 | `registry_confident_domains` | `rules/05_curated.smk` | `registry-confident-domains` | S1, S3, S8 | `results/11_confident_domains.tsv` |
| R15 | `ted_pombe` | `rules/09_external.smk` | `ted-pombe` | S1, S3, S10 | `resources/ted/pombe_ted_segments.tsv` |
| R16 | `ecod_pombe` | `rules/09_external.smk` | `ecod-pombe` | S1, S3, S10 | `resources/ecod/pombe_ecod_segments.tsv` |
| R17 | `dna_bed` | `rules/11_dna.smk` | `dna-bed` | S1, S3, S11 | `data/nucleotide_bed.tsv` |
| R18 | `dna_map` | `rules/11_dna.smk` | `dna-map` | S1, S3, S11 | `12_dna/DNA_domain.tsv` |
| R19 | `dna_organize` | `rules/11_dna.smk` | `dna-organize` | S1, S3, S11 | `12_dna/merged_DNA_domain.tsv` |
| R20 | `dna_genome_bed` | `rules/11_dna.smk` | `dna-genome-bed` | S1, S3, S11 | `12_dna/domains.bed` |

## Appendix B — Script index

| # | Script | Role | Used by |
|---|---|---|---|
| S1 | `scripts/afdb_domains.py` | Executable shell; adds the repo root to `sys.path` and calls `src.cli.main` | R1–R5, R7–R20 |
| S2 | `scripts/compare_afdb_pombase.py` | Standalone AFdb-vs-PomBase comparison; writes tables, report, and workbook | R6 |
| S3 | `src/cli.py` | Typer application; one subcommand per stage, exit-code handling | R1–R5, R7–R20 |
| S4 | `src/files.py` | Config loading, logging, TSV read/write, required-column checks | All CLI commands |
| S5 | `src/registry.py` | Structure selection, sequence matching, registry assembly | R7 |
| S6 | `src/pae.py` | PAE parsing/clustering, `.domain` files, domain table | R8, R9 |
| S7 | `src/plddt.py` | pLDDT extraction, interval utilities, domain statistics | R10, R14 |
| S8 | `src/segments.py` | Domain merge, type annotation, confident-domain caller | R11, R12, R14 |
| S9 | `src/unify.py` | Segment loading, PomBase annotation, unified table | R13 |
| S10 | `src/external.py` | PomBase subsets of the TED and ECOD dumps | R15, R16 |
| S11 | `src/dna.py` | GFF parsing, nucleotide projection, merged BED | R17–R20 |
| S12 | `src/afdb.py` | CIF/PAE scanning and the AFdb entry inventory | R1–R3 |
| S13 | `src/uniprot.py` | UniProt REST fetch and refresh assembly | R4, R5 |

## Appendix C — Dataset declaration index

Declarations live in the `datasets` block of `config/config.yaml` and are
loaded by `src/files.py`. The registry reads all four entries. The
`AFDB_20251107` entry was dropped: nothing referenced it, and the registry
reads v6 through `config/config.yaml`.

| # | Dataset | Source version |
|---|---|---|
| D1 | `AFDB_v6` | `v6` |
| D2 | `AFDB_v4` | `AFDB_v4_20230708` |
| D3 | `LocalAF2_PomBaseDiffAFdb` | `LocalAF2_PomBaseDiffAFdb_20240305` |
| D4 | `AF3_LongProt` | `AF3_LongProt_20241101` |
