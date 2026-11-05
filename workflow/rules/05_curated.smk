# Curated analysis chain.
#
# PAE clustering, pLDDT extraction, domain merging, type annotation, and
# unification. Every rule reads the gene structure registry.

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


rule registry_unify_segments:
    input:
        pae_table=rules.registry_assign_type.output.table,
        ted_table=ancient(f"{RESOURCES}/ted/pombe_ted_segments.tsv"),
        ecod_table=ancient(f"{RESOURCES}/ecod/pombe_ecod_segments.tsv"),
        registry=rules.build_registry.output.registry,
    output:
        SEGMENTS_TABLE,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
        dataset="curated",
        version=CURATED_VERSION,
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-unify-segments --pae-table {input.pae_table} "
        "--ted-table {input.ted_table} --ecod-table {input.ecod_table} "
        "--registry {input.registry} --output {output} --dataset {params.dataset} "
        "--pae-version {params.version} --config {params.config}"
