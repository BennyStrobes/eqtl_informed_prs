#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)


gene_file="$1"
trait_name="$2"
gtex_tissue_file="$3"
coloc_output_dir="$4"
job_number="$5"
total_jobs="$6"

source ~/.bash_profile

module load R/4.0.1

echo $trait_name
echo $job_number


Rscript run_coloc.R $gene_file $trait_name $gtex_tissue_file $coloc_output_dir $job_number $total_jobs

