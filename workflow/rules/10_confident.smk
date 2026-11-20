# Independent 3D-AF-Surfer style confident-domain caller (side branch).

rule registry_confident_domains:
    input:
        plddt_table=rules.registry_extract_plddt.output.table,
        registry=rules.build_registry.output.registry,
    output:
        table=CONFIDENT_TABLE,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-confident-domains --plddt-table {input.plddt_table} "
        "--registry {input.registry} --output {output.table} --config {params.config}"
