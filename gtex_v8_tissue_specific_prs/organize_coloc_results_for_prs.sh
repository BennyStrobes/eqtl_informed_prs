#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)



gtex_tissue_file="$1"
trait_name="$2"
coloc_input_dir="$3"
coloc_output_dir="$4"
chrom_num="$5"

source ~/.bash_profile


python3 organize_coloc_results_for_prs.py $gtex_tissue_file $trait_name $coloc_input_dir $coloc_output_dir $chrom_num