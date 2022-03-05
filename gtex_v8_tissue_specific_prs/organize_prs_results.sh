#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-5:00                         # Runtime in D-HH:MM format
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

coloc_threshold="0.1"
method="coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method




coloc_threshold="0.3"
method="coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method




coloc_threshold="0.5"
method="coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method




coloc_threshold="0.7"
method="coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method




coloc_threshold="0.9"
method="coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="causal_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v1_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method

method="mm_v2_coloc"
python3 organize_prs_results.py $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $coloc_threshold $trait_name $analyzed_ukbb_prs_dir $method
