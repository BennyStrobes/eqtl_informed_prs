#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)



gtex_tissue_file="$1"
trait_name="$2"
cafeh_input_dir="$3"
cafeh_output_dir="$4"
chrom_num="$5"
version="$6"

source ~/.bash_profile


python3 organize_multivariate_cafeh_results_for_prs.py $gtex_tissue_file $trait_name $cafeh_input_dir $cafeh_output_dir $chrom_num $version