# Legacy per-dataset input declarations and PAE clustering:
# PAE JSON matrices -> .domain files -> PAE domain table.
#
# Every declaration below feeds the per-dataset chain (rules 06-11) and is
# disabled by default; the canonical run is driven by the registry instead.
# Each entry pins the input directories, the gene table, and the version label
# recorded in the unified segment table.

DATASETS = {
    # AlphaFold DB release 2023-07-08, UP000002485_284812_SCHPO_v4.
    "AFDB_v4": {
        "enabled": False,
        "source_version": "AFDB_v4_20230708",
        "pae_dir": resolve("reference:PAE/UP000002485_284812_SCHPO_v4"),
        "structure_dir": resolve("reference:UP000002485_284812_SCHPO_v4"),
        "gene_table": resolve("project:resources/pombase/2023-07-01/gene_IDs_names_products.tsv"),
    },
    # Local AlphaFold 2 predictions, 2024-03-05 run (the one behind the legacy vDulab tables).
    "vDulab": {
        "enabled": False,
        "source_version": "LocalAF2_PomBaseDiffAFdb_20240305",
        "pae_dir": resolve("reference:predictions/LocalAF2_PomBaseDiffAFdb_20240305/pae"),
        "structure_dir": resolve("reference:predictions/LocalAF2_PomBaseDiffAFdb_20240305/structure"),
        "gene_table": resolve("project:resources/pombase/2023-07-01/gene_IDs_names_products.tsv"),
    },
    # Custom Dulab2 predictions, 2024-09-04.
    "vDulab2": {
        "enabled": False,
        "source_version": "vDulab2_20240904",
        "pae_dir": resolve("raw:20240904_formated_results_vDulab2_noNaN/PAE_files"),
        "structure_dir": resolve("raw:20240904_formated_results_vDulab2_noNaN/vDulab2"),
        "gene_table": resolve("project:resources/pombase/2023-07-01/gene_IDs_names_products.tsv"),
    },
    # AlphaFold 3 predictions, 2024-11-01.
    "AF3": {
        "enabled": False,
        "source_version": "AF3_LongProt_20241101",
        "pae_dir": resolve("reference:predictions/AF3_LongProt_20241101/pae"),
        "structure_dir": resolve("reference:predictions/AF3_LongProt_20241101/structure"),
        "gene_table": resolve("raw:20241101_AF3_output/AF3_gene_IDs_names_products.tsv"),
    },
    # Dulab2 predictions restricted to proteins present in the DY47073 annotation.
    "DY47073": {
        "enabled": False,
        "source_version": "vDulab2_20240904",
        "pae_dir": resolve("raw:20240904_formated_results_vDulab2_noNaN/PAE_files"),
        "structure_dir": resolve("raw:20240904_formated_results_vDulab2_noNaN/vDulab2"),
        "gene_table": resolve("project:resources/pombase/2023-07-01/gene_IDs_names_products.tsv"),
        "protein_subset": resolve("project:resources/dy47073/protein_subset.txt"),
    },
}

ENABLED_DATASETS = sorted(name for name, entry in DATASETS.items() if entry["enabled"])


rule pae_to_domains:
    input:
        pae_dir=lambda wildcards: DATASETS[wildcards.dataset]["pae_dir"],
    output:
        domains=directory(f"{RESULTS}/{{dataset}}/domains"),
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
        subset=lambda wildcards: subset_flag(DATASETS[wildcards.dataset]),
    threads: 8
    resources:
        mem_mb=8000,
    conda:
        CORE_ENV
    shell:
        "python {params.cli} pae-to-domains --pae-dir {input.pae_dir} "
        "--output-dir {output.domains} --jobs {threads} {params.subset} --config {params.config}"


rule domains_to_table:
    input:
        rules.pae_to_domains.output.domains,
    output:
        table=f"{RESULTS}/{{dataset}}/PAE_domain.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} domains-to-table --domains-dir {input} --output {output.table}"
