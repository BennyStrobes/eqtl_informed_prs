#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-7:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)




gene_file="$1"
trait_name="$2"
gtex_tissue_file="$3"
coloc_output_dir="$4"



source ~/.bash_profile


python3 extract_organized_pph4_mat_from_coloc_results.py $gene_file $trait_name $gtex_tissue_file $coloc_output_dir