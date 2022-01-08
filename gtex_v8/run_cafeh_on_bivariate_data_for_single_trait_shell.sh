#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1



trait_name="$1"
gtex_tissue_file="$2"
processed_bivariate_cafeh_input_dir="$3"
bivariate_cafeh_output_dir="$4"



# Number of parallel jobs per tissue
total_jobs="6"



# Loop through tissues
if false; then
sed 1d $gtex_tissue_file | while read tissue_name tissue_sample_size; do
	echo $tissue_name
	gene_file=$processed_bivariate_cafeh_input_dir$trait_name"_"$tissue_name"_processed_gene_list.txt"
	# Loop through parallel jobs for this tissue
	for job_number in $(seq 0 `expr $total_jobs - "1"`); do
		echo $job_number
		sbatch run_cafeh_on_bivariate_data.sh $gene_file $trait_name $tissue_name $bivariate_cafeh_output_dir $job_number $total_jobs
	done
done 
fi

if false; then
for chrom_num in $(seq 1 22); do
	echo $chrom_num
	python organize_cafeh_results_for_prs.py $gtex_tissue_file $trait_name $processed_bivariate_cafeh_input_dir $bivariate_cafeh_output_dir $chrom_num &
done
fi