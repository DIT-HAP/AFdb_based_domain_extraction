# Unified segment table: PAE + TED + ECOD for registry-selected proteins.

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
