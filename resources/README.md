# resources/

Small, PomBase-only inputs that the pipeline reads directly. Everything here
is tracked in git. The raw TED and ECOD dumps live under the reference root
(`/hugedata/YushengYang/AFdb/TED` and `/hugedata/YushengYang/AFdb/ECOD`).

| Path | Built by | Contents |
|---|---|---|
| `pombase/2026-06-01/` | copied from the PomBase release | Gene table, PeptideStats, peptide FASTA, GFF3, UniProt mapping (canonical) |
| `pombase/2023-07-01/`, `pombase/2025-05-01/` | copied from the PomBase FTP cache | Legacy releases for the v1 regression baseline |
| `uniprot/uniprot_triangulation_2026-06-01.tsv` | one-off UniProt REST query | PomBase/AFdb/current-UniProt relationship for 101 non-identical proteins |
| `uniprot/uniprot_resolution.tsv` | one-off UniProt REST query | Primary/secondary accession resolution for 100 AFdb-only entries |
| `ted/pombe_ted_segments.tsv` | `afdb-domains ted-pombe` | TED domains for PomBase proteins, with CATH family |
| `ecod/pombe_ecod_segments.tsv` | `afdb-domains ecod-pombe` | ECOD domains for PomBase proteins |
| `dy47073/protein_subset.txt` | manual | UniProt IDs of the DY47073 protein set (legacy) |

The TED and ECOD subsets are extracted from the raw dumps by the `ted-pombe`
and `ecod-pombe` rules, so re-running the workflow regenerates them.

The two UniProt tables are pinned snapshots queried on 2026-11-04. They are
deliberately committed instead of fetched at run time so the `relationship`
column stays reproducible; UniProt changes over time.
