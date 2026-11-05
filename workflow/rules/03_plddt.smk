# Per-residue pLDDT extraction from AlphaFold structures.

rule extract_plddt:
    input:
        structure_dir=lambda wildcards: DATASETS[wildcards.dataset]["structure_dir"],
    output:
        table=f"{RESULTS}/{{dataset}}/pLDDT_scores.tsv",
    params:
        cli=str(CLI),
        subset=lambda wildcards: subset_flag(DATASETS[wildcards.dataset]),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} extract-plddt --structure-dir {input.structure_dir} "
        "--output {output.table} {params.subset}"
