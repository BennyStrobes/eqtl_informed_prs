#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB
#SBATCH --nodes=1



gene_file="$1"
trait_name="$2"
tissue_name="$3"
bivariate_cafeh_output_dir="$4"

#conda activate cafeh  # activate environment



python run_cafeh_on_bivariate_data.py $gene_file $trait_name $tissue_name $bivariate_cafeh_output_dir

