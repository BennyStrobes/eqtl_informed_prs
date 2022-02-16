#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)



ukbb_prs_dir="$1"
trait_name="$2"
ukbb_pheno_file1="$3"
ukbb_pheno_file2="$4"
ukbb_pheno_file3="$5"
analyzed_ukbb_prs_dir="$6"

source ~/.bash_profile

echo $trait_name


method="coloc"
coloc_threshold="0.5"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="adaptive_prior_coloc"
coloc_threshold="0.5"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method



method="coloc"
coloc_threshold="0.7"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="adaptive_prior_coloc"
coloc_threshold="0.7"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method



method="coloc"
coloc_threshold="0.9"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="adaptive_prior_coloc"
coloc_threshold="0.9"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method



method="coloc"
coloc_threshold="0.95"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="adaptive_prior_coloc"
coloc_threshold="0.95"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="coloc"
coloc_threshold="0.99"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="adaptive_prior_coloc"
coloc_threshold="0.99"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

