# Gene structure registry.
#
# Selects exactly one structure and one PAE matrix per PomBase coding gene from
# the comparison tables and the AFdb entry inventory.

rule build_registry:
    input:
        merged=rules.compare_afdb_pombase.output.merged,
        afdb_entry_roles=rules.compare_afdb_pombase.output.entry_roles,
        afdb_inventory=rules.afdb_entry_table.output.table,
        triangulation=UNIPROT_TRIANGULATION,
    output:
        registry=REGISTRY_TABLE,
        summary=REGISTRY_SUMMARY,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
        peptide=POMBASE_PEPTIDE,
        afdb_v4=resolve(config["resources"]["registry"]["afdb_v4_dir"]),
        vdulab=resolve(config["resources"]["registry"]["vdulab_structure_dir"]),
        vdulab_fasta=resolve(config["resources"]["registry"]["vdulab_fasta_dir"]),
        af3=resolve(config["resources"]["registry"]["af3_structure_dir"]),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} build-registry --merged-table {input.merged} "
        "--afdb-entry-roles {input.afdb_entry_roles} "
        "--afdb-inventory {input.afdb_inventory} "
        "--uniprot-triangulation {input.triangulation} "
        "--pombase-peptide {params.peptide} "
        "--afdb-v4-dir {params.afdb_v4} --vdulab-structure-dir {params.vdulab} "
        "--vdulab-fasta-dir {params.vdulab_fasta} "
        "--af3-structure-dir {params.af3} "
        "--output {output.registry} --summary-output {output.summary} "
        "--config {params.config}"
