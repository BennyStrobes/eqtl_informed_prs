#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)



gene_file="$1"
trait_name="$2"
tissue_name="$3"
bivariate_cafeh_output_dir="$4"
job_number="$5"
total_jobs="$6"

source ~/.bash_profile

conda activate cafeh  # activate environment



python3 finish_running_cafeh_on_bivariate_data.py $gene_file $trait_name $tissue_name $bivariate_cafeh_output_dir $job_number $total_jobs

