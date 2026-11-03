# downloaded results
# AF_db_dir="/data/c/yangyusheng_optimized/resource/AlphaFold_Dataset/20230708_downloaded"
# AF_db_name="UP000002485_284812_SCHPO_v4"
# gene_info="/data/c/yangyusheng_optimized/resource/Pombase_FTP/release/20230701/gene_IDs_names_products.tsv"
# gene_info="../results/coding_genes_AF_vs_pombase_vs_DY47073_summary_20240912_updated.xlsx"
gene_info="../data/0_raw/20241101_AF3_output/AF3_gene_IDs_names_products.tsv"
# PAE_outputFolder="../data/2_pae_to_domain_igraph_consistent_point5/"

# predicted results using local AF
# AF_db_dir="../data/0_raw/20240305_predicted_results"
# AF_db_name="vDulab"
# PAE_outputFolder="../data/4_pae_to_domain_igraph_consistent_point5_for_those_modified_proteins/"
# mkdir -p $PAE_outputFolder

# predicted results for DY47073
# AF_db_dir="../data/0_raw/20240904_formated_results_vDulab2_noNaN"
# AF_db_name="vDulab2"
# PAE_outputFolder="../data/5_pae_to_domain_igraph_consistent_point5_for_DY47073_proteins/"
# mkdir -p $PAE_outputFolder

# AF3
AF_db_dir="../data/0_raw/20241101_AF3_output"
AF_db_name="UP000002485_284812_SCHPO_v4_and_vAF3"
PAE_outputFolder="../data/7_pae_to_domain_igraph_consistent_point5_for_AF3_all_proteins/"
mkdir -p $PAE_outputFolder

# merged_domains
# AF_db_dir="../data/1_AFdb_vDulab_vDulab2_results_merged/"
# AF_db_name="UP000002485_284812_SCHPO_v4_and_vDulab_and_vDulab2"
# PAE_outputFolder="../data/6_pae_to_domain_igraph_consistent_point5_for_DY47073_all_proteins/"
# mkdir -p $PAE_outputFolder

# pombe_gff="/data/c/yangyusheng_optimized/resource/Pombase_FTP/release/20230701/Schizosaccharomyces_pombe_all_chromosomes.gff3"
# pombe_gff="../references/DY47073_v240116_combined_genome_liftoff_addBraker_240322_nuc_PCgene_4crisprDesign_addchr.gff3"

# run PAE_to_domain
# sh batch_extract_PAE_domain.sh \
#     ${AF_db_dir}"/PAE_files" \
#     $PAE_outputFolder

# python extract_pLDDT_scores.py \
#     -p $AF_db_dir/PDB \
#     -o ../results/${AF_db_name}_pLDDT_scores.tsv
#     # -p $AF_db_dir/$AF_db_name \

# transfer the PAE_to_domain to the table result
# python transfer_residue_to_domain.py \
#     -i $PAE_outputFolder \
#     -o ../results/${AF_db_name}_PAE_domain.tsv

python merge_PAE_domains.py \
    -pae ../results/${AF_db_name}_PAE_domain.tsv \
    -pLDDT ../results/${AF_db_name}_pLDDT_scores.tsv \
    -sysID $gene_info \
    -o ../results/${AF_db_name}_reformated_PAE_domain_based_on_pLDDT.tsv

python assign_domain_type.py \
    -i ../results/${AF_db_name}_reformated_PAE_domain_based_on_pLDDT.tsv \
    -o ../results/${AF_db_name}_annotated_PAE_domain_based_on_pLDDT.tsv

# python extract_confident_domain.py \
#     -pLDDT ../results/${AF_db_name}_pLDDT_scores.tsv \
#     -sysID $gene_info \
#     -o ../results/${AF_db_name}_confident_domains_based_on_pLDDT.tsv

# python create_bed_for_each_nucleotide.py \
#     -g $pombe_gff \
#     -c 56 \
#     -o ../data/3_nucleotide_bed/nucletides_for_each_coding_region.bed

# python map_DNA_with_protein_domain.py \
#     -d ../results/${AF_db_name}_annotated_PAE_domain_based_on_pLDDT.tsv \
#     -n ../data/3_nucleotide_bed/nucletides_for_each_coding_region.bed \
#     -c 56 \
#     -o ../results/${AF_db_name}_DNA_domain.tsv

# python organize_DNA_domains.py \
#     -d ../results/${AF_db_name}_DNA_domain.tsv \
#     -o ../results/${AF_db_name}_merged_DNA_domain.tsv