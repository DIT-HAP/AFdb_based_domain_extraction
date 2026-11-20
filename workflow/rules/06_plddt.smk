# Per-residue pLDDT extraction from the registry-selected structures.

rule registry_extract_plddt:
    input:
        registry=rules.build_registry.output.registry,
    output:
        table=PLDDT_TABLE,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-extract-plddt --registry {input.registry} "
        "--output {output.table} --config {params.config}"
