# Dulab2 predictions restricted to proteins present in the DY47073 annotation.
DATASETS["DY47073"] = {
    "enabled": True,
    "version": "vDulab2",
    "source_version": "vDulab2_20240904",
    "pae_dir": resolve("raw:20240904_formated_results_vDulab2_noNaN/PAE_files"),
    "structure_dir": resolve("raw:20240904_formated_results_vDulab2_noNaN/vDulab2"),
    "gene_table": resolve("project:resources/pombase/2023-07-01/gene_IDs_names_products.tsv"),
    "protein_subset": resolve("project:resources/dy47073/protein_subset.txt"),
}
