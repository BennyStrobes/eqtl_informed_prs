#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=5G                         # Memory total in MiB (for all cores)



ukbb_sample_names_dir="$1"
ukbb_download_data="$2"
ukbb_genotype_data="$3"
cafeh_prs_betas_dir="$4"
ukbb_prs_dir="$5"
chrom="$6"
trait_name="$7"
version="$8"


source ~/.bash_profile

sample_names_keep_file=${ukbb_sample_names_dir}keep.non_british_european_122K.txt

echo $chrom
echo $beta_threshold
echo $version


cafeh_prs_betas_file=$cafeh_prs_betas_dir"cafeh_results_"$trait_name"_"$version"_prs_beta_chrom_"${chrom}".txt"
cafeh_prs_weighted_betas_file=$cafeh_prs_betas_dir"cafeh_results_"$trait_name"_"$version"_prs_weighted_beta_chrom_"${chrom}".txt"
ukbb_prs_file_stem=$ukbb_prs_dir$trait_name"_"$version"_cafeh_prs_chrom_"${chrom}"_"


python3 generate_prs.py ${ukbb_genotype_data}UKB_MAF0.001_v3.${chrom}.bgen ${ukbb_download_data}ukb1404_imp_chr1_v3_withdrawn3.sample ${sample_names_keep_file} $cafeh_prs_betas_file $cafeh_prs_weighted_betas_file $ukbb_prs_file_stem
