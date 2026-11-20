# Registry-selected PAE matrices -> .domain files -> PAE domain table.

rule registry_pae_to_domains:
    input:
        registry=rules.build_registry.output.registry,
    output:
        domains=directory(PAE_DOMAINS),
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    threads: 8
    resources:
        mem_mb=8000,
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-pae-to-domains --registry {input.registry} "
        "--output-dir {output.domains} --jobs {threads} --config {params.config}"


rule registry_domains_to_table:
    input:
        rules.registry_pae_to_domains.output.domains,
    output:
        table=PAE_TABLE,
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-domains-to-table --domains-dir {input} --output {output.table}"
