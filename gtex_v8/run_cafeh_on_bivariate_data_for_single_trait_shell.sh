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





tissue_name="Adipose_Subcutaneous"






gene_file=$processed_bivariate_cafeh_input_dir$trait_name"_"$tissue_name"_processed_gene_list.txt"

sh run_cafeh_on_bivariate_data.sh $gene_file $trait_name $tissue_name $bivariate_cafeh_output_dir