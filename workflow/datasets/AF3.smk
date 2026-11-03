# AlphaFold 3 predictions, 2024-11-01.
DATASETS["AF3"] = {
    "enabled": True,
    "version": "vAF3",
    "source_version": "AF3_20241101",
    "pae_dir": resolve("raw:20241101_AF3_output/PAE_files"),
    "structure_dir": resolve("raw:20241101_AF3_output/PDB"),
    "gene_table": resolve("raw:20241101_AF3_output/AF3_gene_IDs_names_products.tsv"),
}
