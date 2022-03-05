#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=6G                         # Memory total in MiB (for all cores)


gene_file="$1"
trait_name="$2"
gtex_tissue_file="$3"
coloc_output_dir="$4"


causal_version="v1"
python3 run_mm_unweighted_coloc_second_pass.py $gene_file $trait_name $gtex_tissue_file $coloc_output_dir $causal_version

causal_version="v2"
python3 run_mm_unweighted_coloc_second_pass.py $gene_file $trait_name $gtex_tissue_file $coloc_output_dir $causal_version