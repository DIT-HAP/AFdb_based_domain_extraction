# Legacy per-dataset input declarations and PAE clustering:
# PAE JSON matrices -> .domain files -> PAE domain table.
#
# Dataset paths, versions, and flags come from the `datasets` block in
# config.yaml; every dataset is disabled by default and the canonical run is
# driven by the registry instead.

DATASETS = {
    name: {
        "enabled": entry.get("enabled", False),
        "source_version": entry["version"],
        "pae_dir": resolve(entry["pae"]),
        "structure_dir": resolve(entry["structure"]),
        "gene_table": resolve(entry["gene_table"]) if entry.get("gene_table") else None,
        "protein_subset": (
            resolve(entry["protein_subset"]) if entry.get("protein_subset") else None
        ),
    }
    for name, entry in config["datasets"].items()
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
