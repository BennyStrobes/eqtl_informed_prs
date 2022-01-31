#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)



ukbb_prs_dir="$1"
trait_name="$2"
ukbb_pheno_file1="$3"
ukbb_pheno_file2="$4"
ukbb_pheno_file3="$5"
analyzed_ukbb_prs_dir="$6"

source ~/.bash_profile


if false; then
beta_thresh="0.05"
weight_version="unweighted"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $beta_thresh $weight_version $trait_name $analyzed_ukbb_prs_dir
fi
beta_thresh="0.05"
weight_version="weighted"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $beta_thresh $weight_version $trait_name $analyzed_ukbb_prs_dir

if false; then
beta_thresh="0.01"
weight_version="unweighted"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $beta_thresh $weight_version $trait_name $analyzed_ukbb_prs_dir

beta_thresh="0.01"
weight_version="weighted"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $beta_thresh $weight_version $trait_name $analyzed_ukbb_prs_dir


beta_thresh="0.005"
weight_version="unweighted"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $beta_thresh $weight_version $trait_name $analyzed_ukbb_prs_dir

beta_thresh="0.005"
weight_version="weighted"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $beta_thresh $weight_version $trait_name $analyzed_ukbb_prs_dir
fi