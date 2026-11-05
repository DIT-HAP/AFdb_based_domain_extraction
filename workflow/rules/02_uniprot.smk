# UniProt snapshot refresh.
#
# These rules are deliberately not part of `rule all`: the canonical pipeline
# reads the pinned snapshots in resources/uniprot/, so a clean run stays
# offline and reproducible. Run one of these targets explicitly after a
# comparison to produce a dated candidate, then review and promote it.

rule refresh_uniprot_triangulation:
    input:
        pairs=f"{COMPARISON}/05_nonidentical_pairs.tsv",
    output:
        table=f"{UNIPROT_RESULTS}/uniprot_triangulation_refresh.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} refresh-uniprot-triangulation "
        "--nonidentical-pairs {input.pairs} --output {output.table}"


rule refresh_uniprot_resolution:
    input:
        afdb_only=f"{COMPARISON}/06_afdb_only.tsv",
        gene_table=POMBASE_GENE_TABLE,
    output:
        table=f"{UNIPROT_RESULTS}/uniprot_resolution_refresh.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} refresh-uniprot-resolution "
        "--afdb-only {input.afdb_only} --gene-table {input.gene_table} "
        "--output {output.table}"
