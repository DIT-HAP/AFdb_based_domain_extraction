# Gene structure registry.
#
# Selects exactly one structure and one PAE matrix per PomBase coding gene from
# the comparison merged table and the prediction/structure directories.

rule build_registry:
    input:
        merged=rules.compare_afdb_pombase.output.merged,
    output:
        registry=REGISTRY_TABLE,
        summary=REGISTRY_SUMMARY,
        stats=REGISTRY_STATS,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} build-registry --merged-table {input.merged} "
        "--output {output.registry} --summary-output {output.summary} "
        "--stats-output {output.stats} "
        "--config {params.config}"
