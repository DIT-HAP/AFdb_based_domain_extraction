# Custom Dulab predictions, 2024-03-05 release (the one behind the legacy vDulab tables).
DATASETS["vDulab"] = {
    "enabled": True,
    "version": "vDulab",
    "source_version": "vDulab_20240305",
    "pae_dir": resolve("raw:20240305_predicted_results/PAE_files"),
    "structure_dir": resolve("raw:20240305_predicted_results/vDulab"),
    "gene_table": resolve("project:resources/pombase/2023-07-01/gene_IDs_names_products.tsv"),
}
