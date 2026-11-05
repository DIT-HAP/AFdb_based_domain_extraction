# Genomic projection: coding nucleotides -> domain assignment -> merged BED.

rule dna_bed:
    input:
        gff=POMBASE_GFF3,
        peptide_stats=POMBASE_PEPTIDE_STATS,
    output:
        f"{DATA}/nucleotide_bed.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} dna-bed --gff {input.gff} "
        "--peptide-stats {input.peptide_stats} --output {output}"


rule dna_map:
    input:
        nucleotide_bed=f"{DATA}/nucleotide_bed.tsv",
        annotated=rules.assign_type.output.table,
    output:
        table=f"{RESULTS}/{{dataset}}/DNA_domain.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} dna-map --nucleotide-bed {input.nucleotide_bed} "
        "--annotated {input.annotated} --output {output.table}"


rule dna_organize:
    input:
        rules.dna_map.output.table,
    output:
        table=f"{RESULTS}/{{dataset}}/merged_DNA_domain.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} dna-organize --input {input} --output {output.table}"


rule dna_genome_bed:
    input:
        rules.dna_organize.output.table,
    output:
        table=f"{RESULTS}/{{dataset}}/domains.bed",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} dna-genome-bed --input {input} --output {output.table}"
