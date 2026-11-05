# Gene structure registry and the canonical curated run.
#
# The comparison stage writes the detailed AFdb-vs-PomBase tables; the registry
# stage selects one structure and PAE matrix per gene; the curated rules run
# PAE clustering, pLDDT extraction, domain merging, and unification from that
# registry alone.
#
# Output paths carry the result-index number from docs/results_index.md.

COMPARISON = f"{RESULTS}/01_comparison"
REGISTRY_TABLE = f"{RESULTS}/02_gene_structure_registry.tsv"
REGISTRY_SUMMARY = f"{RESULTS}/03_registry_summary.tsv"
PAE_DOMAINS = f"{RESULTS}/04_domains"
PAE_TABLE = f"{RESULTS}/05_PAE_domain.tsv"
PLDDT_TABLE = f"{RESULTS}/06_pLDDT_scores.tsv"
REFORMATED_TABLE = f"{RESULTS}/07_reformated_PAE_domain.tsv"
ANNOTATED_TABLE = f"{RESULTS}/08_annotated_PAE_domain.tsv"
SEGMENTS_TABLE = f"{RESULTS}/09_segments.tsv"
COMPARE_SCRIPT = str(PROJECT_ROOT / "scripts" / "compare_afdb_pombase.py")
POMBASE_2026 = str(PROJECT_ROOT / "resources" / "pombase" / "2026-06-01")
UNIPROT_TRIANGULATION = resolve(config["resources"]["uniprot"]["triangulation"])
UNIPROT_RESOLUTION = resolve(config["resources"]["uniprot"]["resolution"])
REGISTRY_CONF = config["resources"]["registry"]
CURATED_VERSION = "registry_2026-06-01"


rule compare_afdb_pombase:
    input:
        afdb_dir=resolve(REGISTRY_CONF["afdb_dir"]),
        pae_dir=resolve(REGISTRY_CONF["afdb_pae_dir"]),
        gene_table=str(PROJECT_ROOT / "resources" / "pombase" / "2026-06-01" / "gene_IDs_names_products.tsv"),
        peptide=str(PROJECT_ROOT / "resources" / "pombase" / "2026-06-01" / "peptide.fa.gz"),
        triangulation=UNIPROT_TRIANGULATION,
        resolution=UNIPROT_RESOLUTION,
    output:
        merged=f"{COMPARISON}/01_merged_table.tsv",
        entries=f"{COMPARISON}/02_afdb_entries.tsv",
        shared=f"{COMPARISON}/03_shared_comparison.tsv",
        by_sequence=f"{COMPARISON}/04_shared_by_sequence.tsv",
        nonidentical=f"{COMPARISON}/05_nonidentical_pairs.tsv",
        afdb_only=f"{COMPARISON}/06_afdb_only.tsv",
        pombase_only=f"{COMPARISON}/07_pombase_only.tsv",
        counts=f"{COMPARISON}/08_summary_counts.tsv",
        report=f"{COMPARISON}/09_report.md",
        workbook=f"{COMPARISON}/10_afdb_vs_pombase.xlsx",
    params:
        script=COMPARE_SCRIPT,
        pombase_dir=POMBASE_2026,
        release=config["resources"]["pombase"]["coding_release"],
    conda:
        CORE_ENV
    shell:
        "python {params.script} --afdb-dir {input.afdb_dir} --pae-dir {input.pae_dir} "
        "--pombase-dir {params.pombase_dir} --output-dir {COMPARISON} "
        "--pombase-release {params.release} "
        "--uniprot-resolution {input.resolution} "
        "--uniprot-triangulation {input.triangulation}"


rule build_registry:
    input:
        merged=rules.compare_afdb_pombase.output.merged,
        entries=rules.compare_afdb_pombase.output.entries,
        triangulation=UNIPROT_TRIANGULATION,
    output:
        registry=REGISTRY_TABLE,
        summary=REGISTRY_SUMMARY,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
        peptide=str(PROJECT_ROOT / "resources" / "pombase" / "2026-06-01" / "peptide.fa.gz"),
        afdb=resolve(REGISTRY_CONF["afdb_dir"]),
        afdb_v4=resolve(REGISTRY_CONF["afdb_v4_dir"]),
        vdulab=resolve(REGISTRY_CONF["vdulab_structure_dir"]),
        vdulab_fasta=resolve(REGISTRY_CONF["vdulab_fasta_dir"]),
        af3=resolve(REGISTRY_CONF["af3_structure_dir"]),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} build-registry --merged-table {input.merged} "
        "--afdb-entries {input.entries} --uniprot-triangulation {input.triangulation} "
        "--pombase-peptide {params.peptide} --afdb-dir {params.afdb} "
        "--afdb-v4-dir {params.afdb_v4} --vdulab-structure-dir {params.vdulab} "
        "--vdulab-fasta-dir {params.vdulab_fasta} "
        "--af3-structure-dir {params.af3} "
        "--output {output.registry} --summary-output {output.summary} "
        "--config {params.config}"


rule registry_pae_to_domains:
    input:
        registry=rules.build_registry.output.registry,
    output:
        domains=directory(PAE_DOMAINS),
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    threads: 8
    resources:
        mem_mb=8000,
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-pae-to-domains --registry {input.registry} "
        "--output-dir {output.domains} --jobs {threads} --config {params.config}"


rule registry_domains_to_table:
    input:
        rules.registry_pae_to_domains.output.domains,
    output:
        table=PAE_TABLE,
    params:
        cli=str(CLI),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-domains-to-table --domains-dir {input} --output {output.table}"


rule registry_extract_plddt:
    input:
        registry=rules.build_registry.output.registry,
    output:
        table=PLDDT_TABLE,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-extract-plddt --registry {input.registry} "
        "--output {output.table} --config {params.config}"


rule registry_merge_domains:
    input:
        pae_table=rules.registry_domains_to_table.output.table,
        plddt_table=rules.registry_extract_plddt.output.table,
        registry=rules.build_registry.output.registry,
    output:
        table=REFORMATED_TABLE,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} registry-merge-domains --pae-table {input.pae_table} "
        "--plddt-table {input.plddt_table} --registry {input.registry} "
        "--output {output.table} --config {params.config}"


rule registry_assign_type:
    input:
        rules.registry_merge_domains.output.table,
    output:
        table=ANNOTATED_TABLE,
    params:
        cli=str(CLI),
        config=str(CONFIG_PATH),
    conda:
        CORE_ENV
    shell:
        "python {params.cli} assign-type --input {input} --output {output.table} "
        "--config {params.config}"


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
