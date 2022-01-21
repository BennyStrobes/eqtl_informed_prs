#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)



pseudotissue_name="$1"
pseudotissue_expression_dir="$2"
pseudotissue_genotype_dir="$3"
pseudotissue_covariate_dir="$4"
pseudotissue_eqtl_dir="$5"

module load R/3.5.1


cis_distance="500000"



for chrom_num in $(seq 1 22); do 

	# Input files
	matrix_eqtl_covariate_file=${pseudotissue_covariate_dir}${pseudotissue_name}"_covariates.txt"
	matrix_eqtl_expression_file=${pseudotissue_expression_dir}${pseudotissue_name}"_normalized_expression_matrix_eqtl_ready_chr"${chrom_num}".txt"
	matrix_eqtl_gene_loc_file=${pseudotissue_expression_dir}${pseudotissue_name}"_gene_location_matrix_eqtl_ready_chr"${chrom_num}".txt"
	matrix_eqtl_genotype_file=${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_dosage_sample_ordered.txt"
	matrix_eqtl_snp_loc_file=${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_snp_location.txt"
	# Ouput file
	output_file=${pseudotissue_eqtl_dir}${pseudotissue_name}"_chr"${chrom_num}"_matrix_eqtl_results.txt"

	# Run eqtl analysis on this chromosome
	Rscript matrix_eqtl_wrapper.R $matrix_eqtl_genotype_file $matrix_eqtl_snp_loc_file $matrix_eqtl_expression_file $matrix_eqtl_gene_loc_file $matrix_eqtl_covariate_file $cis_distance $output_file

done