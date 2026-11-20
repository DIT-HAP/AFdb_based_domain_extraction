# Merge PAE clusters with pLDDT confidence, then assign domain type.

rule registry_merge_domains:
    input:
        pae_table=rules.registry_domains_to_table.output.table,
        plddt_table=rules.registry_extract_plddt.output.table,
        registry=rules.build_registry.output.registry,
    output:
        table=REFORMATED_TABLE,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-merge-domains --pae-table {input.pae_table} "
        "--plddt-table {input.plddt_table} --registry {input.registry} "
        "--output {output.table} --config {params.config}"


rule registry_assign_type:
    input:
        rules.registry_merge_domains.output.table,
    output:
        table=ANNOTATED_TABLE,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} assign-type --input {input} --output {output.table} "
        "--config {params.config}"
