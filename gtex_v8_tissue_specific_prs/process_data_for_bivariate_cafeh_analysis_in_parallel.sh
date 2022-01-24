#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --mem=20GB
#SBATCH --nodes=1




trait_name="$1"
trait_sumstat_file="$2"
gtex_tissue_file="$3"
genotype_reference_panel_dir="$4"
processed_gtex_associations_dir="$5"
processed_bivariate_cafeh_input_dir="$6"
chrom_num="$7"
sample_size="$8"
cafeh_gene_list_file="$9"

source ~/.bash_profile

python3 process_data_for_bivariate_cafeh_analysis_in_parallel.py $trait_name $trait_sumstat_file $gtex_tissue_file $genotype_reference_panel_dir $processed_gtex_associations_dir $processed_bivariate_cafeh_input_dir $chrom_num $sample_size $cafeh_gene_list_file
