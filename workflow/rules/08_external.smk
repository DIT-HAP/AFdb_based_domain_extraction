# PomBase subsets of the TED and ECOD dumps.

rule ted_pombe:
    input:
        ted_dump=TED_DUMP,
        gene_table=ancient(FEATURES_GENE_TABLE),
    output:
        table=f"{EXTERNAL_RESULTS}/pombe_ted_segments.tsv",
        full=f"{EXTERNAL_RESULTS}/pombe_ted_segments_full_columns.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} ted-pombe --ted-dump {input.ted_dump} "
        "--gene-table {input.gene_table} --output {output.table} "
        "--full-output {output.full}"


rule ecod_pombe:
    input:
        ecod_dump=ECOD_DUMP,
        gene_table=ancient(FEATURES_GENE_TABLE),
    output:
        table=f"{EXTERNAL_RESULTS}/pombe_ecod_segments.tsv",
        full=f"{EXTERNAL_RESULTS}/pombe_ecod_segments_full_columns.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} ecod-pombe --ecod-dump {input.ecod_dump} "
        "--gene-table {input.gene_table} --output {output.table} "
        "--full-output {output.full}"
