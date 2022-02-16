#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)




trait_name="$1"
trait_sumstat_file="$2"
gtex_tissue_file="$3"
genotype_reference_panel_dir="$4"
processed_gtex_associations_dir="$5"
processed_coloc_input_dir="$6"
chrom_num="$7"
sample_size="$8"
cafeh_gene_list_file="$9"

source ~/.bash_profile



python3 process_data_for_coloc_analysis_in_parallel.py $trait_name $trait_sumstat_file $gtex_tissue_file $genotype_reference_panel_dir $processed_gtex_associations_dir $processed_coloc_input_dir $chrom_num $sample_size $cafeh_gene_list_file
