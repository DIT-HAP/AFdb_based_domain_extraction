# resources/

Small, PomBase-only inputs that the pipeline reads directly. Everything here
is tracked in git. The raw TED and ECOD dumps live under the reference root
(`/hugedata/YushengYang/AFdb/TED` and `/hugedata/YushengYang/AFdb/ECOD`).

| Path | Built by | Contents |
|---|---|---|
| `pombase/2026-06-01/` | copied from the PomBase release | Gene table, PeptideStats, peptide FASTA, GFF3, UniProt mapping (canonical) |
| `uniprot/uniprot_triangulation_2026-06-01.tsv` | one-off UniProt REST query; refresh via `refresh_uniprot_triangulation` | PomBase/AFdb/current-UniProt relationship for 101 non-identical proteins |
| `uniprot/uniprot_resolution.tsv` | one-off UniProt REST query; refresh via `refresh_uniprot_resolution` | Primary/secondary accession resolution for 100 AFdb-only entries |
| `dy47073/protein_subset.txt` | manual | UniProt IDs of the DY47073 protein set (legacy) |

The TED and ECOD PomBase subsets are pipeline outputs under
`results/08_external/`; they are no longer committed resources.

The two UniProt tables are pinned snapshots queried on 2026-11-04. They are
deliberately committed instead of fetched at run time so the `relationship`
column stays reproducible; UniProt changes over time.

## Refreshing the UniProt snapshots

The refresh rules are outside `rule all` and require a prior comparison run:

```bash
snakemake refresh_uniprot_triangulation --cores 1 --use-conda
snakemake refresh_uniprot_resolution --cores 1 --use-conda
```

They query UniProt once per accession (1 request/s, retry with backoff) and
write minimal-column candidates to `results/02_uniprot/`. Review the
candidates against the committed tables, then promote one by copying it into
`resources/uniprot/`, updating the path in `config/config.yaml`, and
committing both files.
