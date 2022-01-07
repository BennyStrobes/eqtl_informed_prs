#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1


trait_name="$1"
trait_sumstat_file="$2"
gtex_tissue_file="$3"
gtex_cafeh_data="$4"
processed_gtex_associations_dir="$5"
processed_bivariate_cafeh_input_dir="$6"
sample_size="$7"

if false; then
for chrom_num in $(seq 1 22); do 
	sbatch process_data_for_bivariate_cafeh_analysis_in_parallel.sh $trait_name $trait_sumstat_file $gtex_tissue_file $gtex_cafeh_data $processed_gtex_associations_dir $processed_bivariate_cafeh_input_dir $chrom_num $sample_size
done
fi
python merge_bivariate_cafeh_processed_gene_lists.py $trait_name $gtex_tissue_file $processed_bivariate_cafeh_input_dir
