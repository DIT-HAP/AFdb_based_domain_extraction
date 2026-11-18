# AFdb vs PomBase comparison.
#
# Consumes the AFdb entry inventory and the pinned UniProt snapshots, and
# writes the detailed comparison tables, report, and workbook.

AFDB_SOURCE = resolve(config["resources"]["registry"]["afdb_dir"])


rule compare_afdb_pombase:
    input:
        afdb_entries=rules.afdb_entry_table.output.table,
        gene_table=POMBASE_GENE_TABLE,
        peptide=POMBASE_PEPTIDE,
        triangulation=UNIPROT_TRIANGULATION,
        resolution=UNIPROT_RESOLUTION,
    output:
        merged=f"{COMPARISON}/01_merged_table.tsv",
        entry_roles=f"{COMPARISON}/02_afdb_entry_roles.tsv",
        shared=f"{COMPARISON}/03_shared_comparison.tsv",
        by_sequence=f"{COMPARISON}/04_shared_by_sequence.tsv",
        nonidentical=f"{COMPARISON}/05_nonidentical_pairs.tsv",
        afdb_only=f"{COMPARISON}/06_afdb_only.tsv",
        pombase_only=f"{COMPARISON}/07_pombase_only.tsv",
        counts=f"{COMPARISON}/08_summary_counts.tsv",
        report=f"{COMPARISON}/09_report.md",
        workbook=f"{COMPARISON}/10_afdb_vs_pombase.xlsx",
        coding_genes=f"{COMPARISON}/11_pombase_coding_genes.tsv",
        raw_peptides=f"{COMPARISON}/12_pombase_peptides.tsv",
    params:
        script=str(COMPARE_SCRIPT),
        pombase_dir=POMBASE_2026,
        release=config["resources"]["pombase"]["coding_release"],
        afdb_source=AFDB_SOURCE,
    conda:
        CORE_ENV
    shell:
        "python {params.script} --afdb-entries {input.afdb_entries} "
        "--afdb-source {params.afdb_source} "
        "--pombase-dir {params.pombase_dir} --output-dir {COMPARISON} "
        "--pombase-release {params.release} "
        "--uniprot-resolution {input.resolution} "
        "--uniprot-triangulation {input.triangulation}"
