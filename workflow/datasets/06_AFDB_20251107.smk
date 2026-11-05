# AlphaFold DB release 2025-11-07, UP000002485_284812_SCHPO_v6.
# Configured for the next upgrade; not part of the current acceptance run.
DATASETS["AFDB_20251107"] = {
    "enabled": False,
    "version": "v6",
    "source_version": "AFDB_v6_20251107",
    "pae_dir": resolve("cache:AlphaFold_Dataset/20251107_downloaded/PAE_files"),
    "structure_dir": resolve("cache:AlphaFold_Dataset/20251107_downloaded/UP000002485_284812_SCHPO_v6"),
    "gene_table": resolve("project:resources/pombase/2025-05-01/gene_IDs_names_products.tsv"),
}
