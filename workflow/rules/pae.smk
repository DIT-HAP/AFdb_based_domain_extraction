# PAE clustering: JSON matrices -> .domain files -> PAE domain table.

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
