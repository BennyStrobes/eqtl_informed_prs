#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)

source ~/.bash_profile



pseudotissue_info_file="$1"
pseudotissue_sample_names_dir="$2"
gtex_tpm_expression_dir="$3"
gene_annotation_file="$4"
pseudotissue_expression_dir="$5"

# Loop through pseudotissues in pseudotissue info file
sed 1d $pseudotissue_info_file | while read pseudotissue sample_size sample_repeat composit_tissue_string; do

	echo ${pseudotissue}"\t"${composit_tissue_string}
	python3 generate_pseudotissue_expression_matrix.py $pseudotissue $composit_tissue_string $pseudotissue_sample_names_dir $gtex_tpm_expression_dir $gene_annotation_file $pseudotissue_expression_dir
	expression_file=${pseudotissue_expression_dir}${pseudotissue}"_normalized_expression_autosomes_protein_coding_and_linc.txt"


	# Required for matrix eqtl ready
	## Put in required format
	## Seperate into different chromosomes
	matrix_eqtl_expression_output_root=${pseudotissue_expression_dir}${pseudotissue}"_normalized_expression_matrix_eqtl_ready_"
	matrix_eqtl_gene_location_output_root=${pseudotissue_expression_dir}${pseudotissue}"_gene_location_matrix_eqtl_ready_"
	python3 make_expression_data_matrix_eqtl_ready.py $expression_file $matrix_eqtl_expression_output_root $matrix_eqtl_gene_location_output_root
done




