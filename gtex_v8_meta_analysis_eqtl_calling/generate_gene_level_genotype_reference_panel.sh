#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=60G                         # Memory total in MiB (for all cores)

source ~/.bash_profile

chrom_num="$1"
cross_tissue_gene_list_file="$2"
gtex_genotype_dir="$3"
genotype_reference_panel_dir="$4"

plink --bfile ${gtex_genotype_dir}"GTEx_v8_genotype_EUR."${chrom_num} --geno 0.0 --maf .025 --recode vcf --out ${genotype_reference_panel_dir}"GTEx_v8_genotype_EUR_"${chrom_num}


python3 convert_vcf_file_to_dosage_file.py ${genotype_reference_panel_dir}"GTEx_v8_genotype_EUR_"${chrom_num}


rm ${genotype_reference_panel_dir}"GTEx_v8_genotype_EUR_"${chrom_num}".vcf"
rm ${genotype_reference_panel_dir}"GTEx_v8_genotype_EUR_"${chrom_num}".nosex"
rm ${genotype_reference_panel_dir}"GTEx_v8_genotype_EUR_"${chrom_num}".log"


genotype_dosage_file=${genotype_reference_panel_dir}"GTEx_v8_genotype_EUR_"${chrom_num}"_dosage.txt"
python3 extract_gene_level_genotype_reference_panels_on_single_chromosome.py $chrom_num $genotype_dosage_file $cross_tissue_gene_list_file $genotype_reference_panel_dir
