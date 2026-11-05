# Unified segment table per dataset and across datasets.

rule unify_segments:
    input:
        pae_table=rules.assign_type.output.table,
        ted_table=f"{RESOURCES}/ted/pombe_ted_segments.tsv",
        ecod_table=f"{RESOURCES}/ecod/pombe_ecod_segments.tsv",
        gene_table=FEATURES_GENE_TABLE,
        stats_table=FEATURES_PEPTIDE_STATS,
    output:
        f"{RESULTS}/{{dataset}}/segments.tsv",
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
        dataset=lambda wildcards: wildcards.dataset,
        version=lambda wildcards: DATASETS[wildcards.dataset]["source_version"],
    wildcard_constraints:
        dataset="(?!combined$)[A-Za-z0-9_.-]+",
    conda:
        CORE_ENV
    shell:
        "python {params.cli} unify-segments --pae-table {input.pae_table} "
        "--ted-table {input.ted_table} --ecod-table {input.ecod_table} "
        "--gene-table {input.gene_table} --stats-table {input.stats_table} "
        "--output {output} --dataset {params.dataset} --pae-version {params.version} "
        "--config {params.config}"


rule combine_segments:
    input:
        tables=expand(f"{RESULTS}/{{dataset}}/segments.tsv", dataset=ENABLED_DATASETS),
    output:
        f"{RESULTS}/combined/segments.tsv",
    params:
        cli=str(CLI),
        flags=lambda wildcards, input: " ".join(f"--input {path}" for path in input.tables),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} combine-segments {params.flags} --output {output}"
