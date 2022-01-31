#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-18:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)


gene_file="$1"
trait_name="$2"
gtex_tissue_file="$3"
cafeh_output_dir="$4"
job_number="$5"
total_jobs="$6"
version="$7"

source ~/.bash_profile

conda activate cafeh  # activate environment


echo $trait_name
echo $job_number



python3 run_cafeh_on_multivariate_data.py $gene_file $trait_name $gtex_tissue_file $cafeh_output_dir $job_number $total_jobs $version

