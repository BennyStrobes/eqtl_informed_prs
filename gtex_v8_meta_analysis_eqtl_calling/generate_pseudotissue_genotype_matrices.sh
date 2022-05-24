#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)

source ~/.bash_profile

pseudotissue_name="$1"
pseudotissue_individual_names="$2"
tissue_sample_names="$3"
gtex_genotype_dir="$4"
pseudotissue_genotype_dir="$5"

echo $pseudotissue_name

for chrom_num in $(seq 1 22); do 
	echo $chrom_num
	# Filter to individuals in this tissue (and re-order columns so they line up with order of expression data)
	plink --bfile ${gtex_genotype_dir}"GTEx_v8_genotype_EUR_1kg_overlap_maf_05_"${chrom_num} --keep-allele-order --keep ${pseudotissue_individual_names} --indiv-sort f ${pseudotissue_individual_names} --make-bed --out ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}
	
	# Get Allele frequency of the snps in this tissue
	plink --bfile ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num} --keep-allele-order --freq --out ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}
	# Print what the minimum allele frequency is (more just double checking we are not getting really low af snps using this approach)
	python3 print_min_af_from_plink_frq_file.py ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}".frq"

	# Convert to VCF in order to convert to matrix eqtl genotype data
	plink --bfile ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num} --keep-allele-order --recode vcf --out ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}
	python3 convert_vcf_file_to_dosage_file.py ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}

	# Get snp location file for matrix eqtl
	genotype_file=${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_dosage.txt"
	snp_loc_file=${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}"_snp_location.txt"
	python3 generate_snp_location_file_from_genotype_file.py $genotype_file $snp_loc_file

	# Remove unecessary files
	rm ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}".log"
	rm ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}".nosex"
	rm ${pseudotissue_genotype_dir}${pseudotissue_name}"_GTEx_v8_genotype_EUR_"${chrom_num}".vcf"	
done



