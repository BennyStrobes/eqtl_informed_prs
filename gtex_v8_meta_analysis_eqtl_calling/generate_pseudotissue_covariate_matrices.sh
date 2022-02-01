#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                         # Memory total in MiB (for all cores)

source ~/.bash_profile


tissue_info_file="$1"
pseudotissue_sample_names_dir="$2"
gtex_covariate_dir="$3"
pseudotissue_covariate_dir="$4"


# Loop through pseudotissues in pseudotissue info file
sed 1d $tissue_info_file | while read tissue_name sample_size pseudotissue_name; do

	echo ${tissue_name}
	pseudotissue_sample_names=${pseudotissue_sample_names_dir}${tissue_name}"_individual_names.txt"

	pseudotissue_covariate_output_file=${pseudotissue_covariate_dir}${tissue_name}"_covariates.txt"
	python3 generate_pseudotissue_covariate_matrix.py $tissue_name $pseudotissue_sample_names $gtex_covariate_dir $pseudotissue_covariate_output_file

done
