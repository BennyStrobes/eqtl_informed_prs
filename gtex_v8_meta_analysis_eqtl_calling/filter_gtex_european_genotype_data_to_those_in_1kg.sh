#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-10:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=60G                         # Memory total in MiB (for all cores)

source ~/.bash_profile




gtex_genotype_dir="$1"
ref_1kg_genotype_dir="$2"
gtex_genotype_1kg_overlap_dir="$3"



for chrom_num in $(seq 1 22); do 
	echo $chrom_num
	ref_1kg_gtex_formatted_bim=$gtex_genotype_1kg_overlap_dir"ref_1kg_chrom_"${chrom_num}"_snps_gtex_formatted.bim"
	python3 convert_1kg_bim_file_to_gtex_snp_id_format.py $ref_1kg_genotype_dir"1000G.EUR.hg38."${chrom_num}".bim" $ref_1kg_gtex_formatted_bim

	echo ${gtex_genotype_dir}"GTEx_v8_genotype_EUR."${chrom_num}
	plink --bfile ${gtex_genotype_dir}"GTEx_v8_genotype_EUR."${chrom_num} --geno 0.0 --maf .05 --extract $ref_1kg_gtex_formatted_bim --make-bed --out ${gtex_genotype_1kg_overlap_dir}"GTEx_v8_genotype_EUR_1kg_overlap_maf_05_"${chrom_num}
done

