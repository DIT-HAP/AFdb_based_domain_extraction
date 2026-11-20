# Legacy output manifest

Inventory of everything in the v1 `results/` tree and its status in v2. For the
numbered inventory of the active outputs, see [results_index.md](results_index.md).

The v1 repository is tagged `v1-legacy` and remains the read-only reference
baseline. Nothing is deleted from it.

## Canonical equivalents

The legacy per-dataset branch was removed on 2026-11-19, so the v1 per-dataset
tables are no longer reproduced directly. The registry chain produces one
canonical table per stage; v1 remains available in git history.

| Group | v1 files | Canonical producer | Canonical output |
|---|---|---|---|
| PAE domain tables | `<ds>_PAE_domain.tsv` | `registry-domains-to-table` | `results/06_PAE_domain.tsv` |
| pLDDT tables | `<ds>_pLDDT_scores.tsv` | `registry-extract-plddt` | `results/07_pLDDT_scores.tsv` |
| Reformatted domains | `<ds>_reformated_PAE_domain_based_on_pLDDT.tsv` | `registry-merge-domains` | `results/08_reformated_PAE_domain.tsv` |
| Annotated domains | `<ds>_annotated_PAE_domain_based_on_pLDDT.tsv` | `assign-type` | `results/09_annotated_PAE_domain.tsv` |
| Confident domains | `<ds>_confident_domains_based_on_pLDDT.tsv` | `registry-confident-domains` | `results/11_confident_domains.tsv` |
| TED PomBase subset | `references/TED/pombe_ted_domains_detailed.tsv` | `ted-pombe` | `resources/ted/pombe_ted_segments.tsv` |
| ECOD PomBase subset | `references/ECOD/pombe_ECOD_domains.tsv` | `ecod-pombe` | `resources/ecod/pombe_ecod_segments.tsv` |

The v1 datasets covered were `AFDB_v4`, `LocalAF2_PomBaseDiffAFdb`,
`LocalAF2_DY47073Annotation`, `AF3_LongProt`, `DY47073_ProteinSubset`.

## New in v2

| Group | v2 files | Producer |
|---|---|---|
| AFdb structure inventory | `results/01_afdb/*` | `afdb-*` |
| AFdb/PomBase comparison | `results/03_comparison/*` | `compare-afdb-pombase` |
| Gene structure registry | `results/04_registry/01_gene_structure_registry.tsv`, `results/04_registry/03_registry_stats.tsv` | `build-registry` |
| Curated unified segments | `results/10_segments.tsv` | `registry-unify-segments` |
| Confident domains | `results/11_confident_domains.tsv` | `registry-confident-domains` |
| Resource provenance | `resources/manifest.tsv` | `fetch-resources` |
| DNA projections | `results/12_dna/DNA_domain.tsv`, `merged_DNA_domain.tsv`, `domains.bed` | `dna-*` |

## Archived, not reproduced

| Group | v1 files | Reason |
|---|---|---|
| Coding-gene summaries | `coding_genes_AF_vs_pombase_summary*.xlsx` (5 dated variants), `compare_AF_pombase.csv` | Exploratory notebook output; superseded by the segment tables |
| AFdb vs vDulab comparison | `AFdb_vDulab_domain_comparison.xlsx` | Notebook output; the unified table carries the same information |
| Confident-domain comparison | `compare_domain_results.xlsx` | Notebook output |
| pLDDT snapshots | `pLDDTs.xlsx`, `20240904_updated_pLDDTs.xlsx`, `20240912_updated_pLDDTs.xlsx`, `20240929_updated_pLDDTs.xlsx` | Superseded by `<ds>_pLDDT_scores.tsv` |
| DY47073 analysis | `duplicated_proteins_in_DY47073.xlsx`, `coding_genes_AF_vs_pombase_vs_DY47073_summary*.xlsx`, `results/Fasta_for_customed_prediction/`, `results/DY47073_fasta_for_customed_prediction/` | Notebook analysis; prediction FASTA sets are downstream inputs, not pipeline outputs |
| DIT_HAP / CKHS figures | `results/domain_comparison/*.pdf` (16 files) | Out of scope by decision (2026-11-03); data and analysis removed |

## Known differences from v1

1. **Duplicate legacy rows.** v1 contains duplicated rows for three proteins
   (`O14460`, `P09322`, `P09988`) in the confident-domain table and in the
   AFDB v4 reformatted/annotated tables (9 extra rows). v2 de-duplicates; the
   residue intervals and `fraction_gt_70` match exactly.
2. **Rounding drift.** The confident-domain table drifts by ±0.01
   `mean_pLDDT` on 45 rows. `LocalAF2_PomBaseDiffAFdb` and
   `LocalAF2_DY47073Annotation` drift on 2–3 rows each, and AFDB v4 drifts on
   68 of 17,851 rows. All residue intervals, lengths, `fraction_gt_70`,
   locations and domain types match; AF3_LongProt matches exactly on every
   field.
3. **Duplicate legacy inputs.** The v1 `vDulab` tables (now
   LocalAF2_PomBaseDiffAFdb) came from `20240305_predicted_results`, not from
   the later renamed release. v2 pins the 2024-03-05 directory; the 45
   sequence-consistent models are copied to
   `predictions/LocalAF2_PomBaseDiffAFdb_20240305` on hugedata.
4. **Column order.** v2 writes `domain_id/domain_residues` before
   `protein_length`; values are identical.
5. **UniProt fallback in `Systematic ID`.** The v1 `merge_PAE_domains.py`
   filled unmapped proteins with their UniProt accession. The legacy branch
   was removed on 2026-11-19; the registry tables join on the registry's own
   gene ids, so the fallback no longer exists.
