# PomBase subsets of the TED and ECOD dumps.

rule ted_pombe:
    input:
        ted_dump=TED_DUMP,
        gene_table=ancient(FEATURES_GENE_TABLE),
    output:
        f"{RESOURCES}/ted/pombe_ted_segments.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} ted-pombe --ted-dump {input.ted_dump} "
        "--gene-table {input.gene_table} --output {output}"


rule ecod_pombe:
    input:
        ecod_dump=ECOD_DUMP,
        gene_table=ancient(FEATURES_GENE_TABLE),
    output:
        f"{RESOURCES}/ecod/pombe_ecod_segments.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} ecod-pombe --ecod-dump {input.ecod_dump} "
        "--gene-table {input.gene_table} --output {output}"
