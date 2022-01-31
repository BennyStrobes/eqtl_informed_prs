#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1



trait_name="$1"
gtex_tissue_file="$2"
cafeh_input_dir="$3"
cafeh_output_dir="$4"
version="$5"

source ~/.bash_profile

# Number of parallel jobs
total_jobs="20"

gene_file=$cafeh_input_dir$trait_name"_processed_gene_list.txt"



# Loop through parallel jobs for this tissue
if false; then
for job_number in $(seq 0 `expr $total_jobs - "1"`); do
	sbatch run_cafeh_on_multivariate_data.sh $gene_file $trait_name $gtex_tissue_file $cafeh_output_dir $job_number $total_jobs $version
done
fi


if false; then
for chrom_num in $(seq 1 22); do
	echo $chrom_num
	sbatch organize_multivariate_cafeh_results_for_prs.sh $gtex_tissue_file $trait_name $cafeh_input_dir $cafeh_output_dir $chrom_num $version
done
fi

if false; then
python3 print_number_of_cafeh_components_per_tissue.py $gtex_tissue_file $trait_name $cafeh_output_dir $version
fi
