#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 4-12:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=6G                         # Memory total in MiB (for all cores)

pseudotissue_name="$1"
pseudotissue_expression_dir="$2"
pseudotissue_genotype_dir="$3"
pseudotissue_covariate_dir="$4"
pseudotissue_sample_names_dir="$5"
pseudotissue_eqtl_dir="$6"
chrom_num="$7"

source ~/.bash_profile

module load R/3.6.1

cis_distance="500000"

echo ${pseudotissue_name}"_"${chrom_num}
	# Input files
	matrix_eqtl_covariate_file=${pseudotissue_covariate_dir}${pseudotissue_name}"_covariates.txt"
	matrix_eqtl_expression_file=${pseudotissue_expression_dir}${pseudotissue_name}"_normalized_expression_matrix_eqtl_ready_chr"${chrom_num}".txt"
	matrix_eqtl_gene_loc_file=${pseudotissue_expression_dir}${pseudotissue_name}"_gene_location_matrix_eqtl_ready_chr"${chrom_num}".txt"
	matrix_eqtl_genotype_file=${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_dosage_sample_ordered.txt"
	matrix_eqtl_snp_loc_file=${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_snp_location.txt"
	matrix_eqtl_sample_repeat_file=${pseudotissue_sample_names_dir}${pseudotissue_name}"_sample_repeate_info.txt"
	# Ouput file
	output_file=${pseudotissue_eqtl_dir}${pseudotissue_name}"_chr"${chrom_num}"_lmm_eqtl_results.txt"


	Rscript lmm_eqtl.R $matrix_eqtl_genotype_file $matrix_eqtl_snp_loc_file $matrix_eqtl_expression_file $matrix_eqtl_gene_loc_file $matrix_eqtl_covariate_file $matrix_eqtl_sample_repeat_file $cis_distance $output_file