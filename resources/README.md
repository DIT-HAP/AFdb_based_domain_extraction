# resources/

Small, PomBase-only inputs that the pipeline reads directly. Everything here
is tracked in git. The raw TED and ECOD dumps live under the reference root
(`/hugedata/YushengYang/AFdb/TED` and `/hugedata/YushengYang/AFdb/ECOD`).

| Path | Built by | Contents |
|---|---|---|
| `pombase/2023-07-01/` | copied from the PomBase FTP cache | Gene table, PeptideStats, GFF3 for the coding release |
| `pombase/2025-05-01/` | copied from the PomBase FTP cache | Gene table and PeptideStats for TED/ECOD annotation |
| `ted/pombe_ted_segments.tsv` | `afdb-domains ted-pombe` | TED domains for PomBase proteins, with CATH family |
| `ecod/pombe_ecod_segments.tsv` | `afdb-domains ecod-pombe` | ECOD domains for PomBase proteins |
| `dy47073/protein_subset.txt` | manual | UniProt IDs of the DY47073 protein set |

The TED and ECOD subsets are extracted from the raw dumps by the `ted-pombe`
and `ecod-pombe` rules, so re-running the workflow regenerates them.
