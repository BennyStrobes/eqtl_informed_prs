#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)

source ~/.bash_profile


downsampled_tissue_info_file="$1"
gtex_downsampled_individuals_dir="$2"
pseudotissue_sample_names_dir="$3"



python3 generate_pseudotissue_sample_names.py $downsampled_tissue_info_file $gtex_downsampled_individuals_dir $pseudotissue_sample_names_dir
