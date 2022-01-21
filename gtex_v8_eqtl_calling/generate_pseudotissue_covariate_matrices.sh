#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                         # Memory total in MiB (for all cores)

source ~/.bash_profile


pseudotissue_info_file="$1"
pseudotissue_sample_names_dir="$2"
pseudotissue_expression_dir="$3"
gtex_covariate_dir="$4"
pseudotissue_covariate_dir="$5"




# Loop through pseudotissues in pseudotissue info file
sed 1d $pseudotissue_info_file | while read pseudotissue sample_size sample_repeat composit_tissue_string; do

	echo ${pseudotissue}"\t"${composit_tissue_string}


	# Input files
	pseudotissue_sample_names=${pseudotissue_sample_names_dir}${pseudotissue}"_sample_names.txt"
	pseudotissue_expr_pc_file=${pseudotissue_expression_dir}${pseudotissue}"_normalized_expression_pcs.txt"
	# Output file
	pseudotissue_covariate_output_file=${pseudotissue_covariate_dir}${pseudotissue}"_covariates.txt"


	python3 generate_pseudotissue_covariate_matrix.py $pseudotissue $composit_tissue_string $pseudotissue_sample_names $pseudotissue_expr_pc_file $gtex_covariate_dir $pseudotissue_covariate_output_file



done
