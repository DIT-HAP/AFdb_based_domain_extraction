# Protein domain merging, annotation, and the independent confident caller.

rule merge_domains:
    input:
        pae_table=rules.domains_to_table.output.table,
        plddt_table=rules.extract_plddt.output.table,
        gene_table=lambda wildcards: DATASETS[wildcards.dataset]["gene_table"],
    output:
        table=f"{RESULTS}/{{dataset}}/reformated_PAE_domain.tsv",
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} merge-domains --pae-table {input.pae_table} "
        "--plddt-table {input.plddt_table} --gene-table {input.gene_table} "
        "--output {output.table} --config {params.config}"


rule assign_type:
    input:
        rules.merge_domains.output.table,
    output:
        table=f"{RESULTS}/{{dataset}}/annotated_PAE_domain.tsv",
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} assign-type --input {input} --output {output.table} --config {params.config}"


rule confident_domains:
    input:
        plddt_table=rules.extract_plddt.output.table,
        gene_table=lambda wildcards: DATASETS[wildcards.dataset]["gene_table"],
    output:
        table=f"{RESULTS}/{{dataset}}/confident_domains.tsv",
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} confident-domains --plddt-table {input.plddt_table} "
        "--gene-table {input.gene_table} --output {output.table} --config {params.config}"
