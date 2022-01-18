#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)


ukbb_sample_names_dir="$1"
ukbb_download_data="$2"
cafeh_prs_betas_dir="$3"
ukbb_prs_dir="$4"
chrom="$5"
beta_threshold="$6"
trait_name="$7"


source ~/.bash_profile

sample_names_keep_file=${ukbb_sample_names_dir}keep.non_british_european_122K.txt

echo $chrom
echo $beta_threshold


cafeh_prs_betas_file=$cafeh_prs_betas_dir"cafeh_results_"$trait_name"_prs_beta_chrom_"${chrom}".txt"
ukbb_prs_file=$ukbb_prs_dir$trait_name"_cafeh_prs_beta_threshold_"${beta_threshold}"_chrom_"${chrom}".txt"

python3 generate_prs.py ${ukbb_download_data}ukb_imp_chr${chrom}_v3.bgen ${ukbb_download_data}ukb1404_imp_chr1_v3_withdrawn3.sample ${sample_names_keep_file} $cafeh_prs_betas_file $ukbb_prs_file $beta_threshold
