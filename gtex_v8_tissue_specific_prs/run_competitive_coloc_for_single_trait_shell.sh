#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1



trait_name="$1"
gtex_tissue_file="$2"
coloc_input_dir="$3"
coloc_output_dir="$4"


# Number of parallel jobs
total_jobs="20"

gene_file=$coloc_input_dir$trait_name"_processed_gene_list.txt"


if false; then
for job_number in $(seq 0 `expr $total_jobs - "1"`); do
	sbatch run_competitive_coloc.sh $gene_file $trait_name $gtex_tissue_file $coloc_output_dir $job_number $total_jobs
done
fi


if false; then
for chrom_num in $(seq 1 22); do
	echo $chrom_num
	sbatch organize_competitive_coloc_results_for_prs.sh $gtex_tissue_file $trait_name $coloc_input_dir $coloc_output_dir $chrom_num
done
fi

if false; then
python3 print_number_of_coloc_components_per_tissue.py $gtex_tissue_file $trait_name $coloc_output_dir
fi
