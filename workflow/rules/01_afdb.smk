# AFdb structure inventory.
#
# Scans the pinned AFdb v6 CIF headers and PAE filenames once, then joins them
# into the single entry table consumed by the comparison and registry stages.

AFDB_CIF_DIR = resolve(config["datasets"]["AFDB_v6"]["structure"])
AFDB_PAE_DIR = resolve(config["datasets"]["AFDB_v6"]["pae"])


rule afdb_cif_entries:
    input:
        afdb_dir=AFDB_CIF_DIR,
    output:
        table=f"{AFDB_RESULTS}/01_afdb_cif_entries.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} afdb-cif-entries --afdb-dir {input.afdb_dir} "
        "--output {output.table}"


rule afdb_pae_index:
    input:
        pae_dir=AFDB_PAE_DIR,
    output:
        table=f"{AFDB_RESULTS}/02_afdb_pae_index.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} afdb-pae-index --pae-dir {input.pae_dir} "
        "--output {output.table}"


rule afdb_entry_table:
    input:
        cif_entries=rules.afdb_cif_entries.output.table,
        pae_index=rules.afdb_pae_index.output.table,
    output:
        table=f"{AFDB_RESULTS}/03_afdb_entries.tsv",
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} afdb-entry-table --cif-entries {input.cif_entries} "
        "--pae-index {input.pae_index} --output {output.table}"
