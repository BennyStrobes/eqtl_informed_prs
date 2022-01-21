#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)

source ~/.bash_profile

pseudotissue_name="$1"
pseudotissue_individual_names="$2"
pseudotime_sample_names="$3"
gtex_genotype_dir="$4"
pseudotissue_genotype_dir="$5"

echo $pseudotissue_name



for chrom_num in $(seq 1 22); do 

	echo $chrom_num
	plink --bfile ${gtex_genotype_dir}"GTEx_v8_genotype_EUR."${chrom_num} --keep ${pseudotissue_individual_names} --geno 0.0 --maf .01 --recode vcf --out ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}



	python3 convert_vcf_file_to_dosage_file.py ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}


	python3 reorder_vcf_dosage_columns.py ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num} $pseudotissue_individual_names $pseudotime_sample_names
	genotype_file=${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_dosage_sample_ordered.txt"

	snp_loc_file=${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_snp_location.txt"
	python3 generate_snp_location_file_from_genotype_file.py $genotype_file $snp_loc_file


	rm ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_dosage.txt"
	rm ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}".log"
	rm ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}".nosex"
	rm ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}".vcf"
done
