#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)


ukbb_sample_names_dir="$1"
ukbb_download_data="$2"
ukbb_genotype_data="$3"
cafeh_prs_betas_dir="$4"
ukbb_prs_dir="$5"
chrom="$6"
trait_name="$7"


source ~/.bash_profile

sample_names_keep_file=${ukbb_sample_names_dir}keep.non_british_european_122K.txt

echo $trait_name
echo $chrom


python3 generate_prs.py ${ukbb_genotype_data}"UKB_MAF0.001_v3."${chrom}".bgen" ${ukbb_download_data}"ukb1404_imp_chr1_v3_withdrawn3.sample" ${sample_names_keep_file} $cafeh_prs_betas_dir $ukbb_prs_dir $chrom $trait_name