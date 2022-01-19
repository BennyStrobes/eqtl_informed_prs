#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=4G                         # Memory total in MiB (for all cores)



ukbb_prs_dir="$1"
ukbb_pheno_file1="$2"
ukbb_pheno_file2="$3"
ukbb_pheno_file3="$4"

source ~/.bash_profile
if false; then
thresh=".005"
echo $thresh
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $thresh


thresh=".05"
echo $thresh
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $thresh

thresh=".1"
echo $thresh
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $thresh


thresh=".5"
echo $thresh
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $thresh

thresh="1"
echo $thresh
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $thresh
fi


thresh=".05"
echo $thresh
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $thresh



