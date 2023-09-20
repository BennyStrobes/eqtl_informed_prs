#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=60G                         # Memory total in MiB (for all cores)

source ~/.bash_profile




gtex_genotype_dir="$1"
gtex_processed_genotype_dir="$2"



for chrom_num in $(seq 1 22); do 
	plink --bfile ${gtex_genotype_dir}"GTEx_v8_genotype_EUR."${chrom_num} --geno 0.0 --maf .05 --keep-allele-order --make-bed --out ${gtex_processed_genotype_dir}"GTEx_v8_genotype_EUR_maf_05_"${chrom_num}
done
