#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-2:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                         # Memory total in MiB (for all cores)

##################
# Input data
##################
# Directory containing GTEx genotype data
# Shared by Tiffany. Created by Huwenbo
gtex_genotype_dir="/n/groups/price/huwenbo/DATA/GTEx_v8/GTEx_v8_genotypes_EUR/"

# Directory containing GTEx covariate data
# Shared by Tiffany. Created by Huwenbo
gtex_covariate_dir="/n/groups/price/huwenbo/DATA/GTEx_v8/GTEx_Analysis_v8_eQTL_covariates/"

# Directory containing GTEx TPM expression
gtex_tpm_expression_dir="/n/groups/price/GTEX/GTEX_V8/TPM/"

# Directory continaing file for each tissue (Downsampled_Individuals_320_{$tissue_name}.txt) containing down-sampled individuals
# Directory created by Tiffany.
gtex_downsampled_individuals_dir="/n/groups/price/tiffany/subpheno/AllGTExTissues_restore/Downsampled_Ind/"

# File created by Tiffany containing info on downsampled tissues
# Shared by Tiffany on slack on Jan 13, 2022
downsampled_tissue_info_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/input_data/TissueGroups_v2.txt"

# GTEx gencode gene annotation file
# Downloaded from https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf on Jan 19 2022
gene_annotation_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/input_data/gencode.v26.GRCh38.genes.gtf"


##################
# Output data
##################
#Root ouptut directory
output_root="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/"
# Directory containing file for each tissue with ordered sample names in pseudotissues
pseudotissue_sample_names_dir=${output_root}"pseudotissue_sample_names/"
# Directory containing expression file for each pseudotissue
pseudotissue_expression_dir=${output_root}"pseudotissue_expression/"
# Directory containing genotype daata
pseudotissue_genotype_dir=${output_root}"pseudotissue_genotype/"
# Directory containing covariate daata
pseudotissue_covariate_dir=${output_root}"pseudotissue_covariates/"
# Directory containing eQTL results
pseudotissue_eqtl_dir=${output_root}"pseudotissue_eqtl_summary_stats/"
# Directory containing gene level genotype reference panel
genotype_reference_panel_dir=${output_root}"gene_level_genotype_reference_panel/"

##################
# Run analysis
##################


###################
# create ordered list of samples for each of the pseudotissues
if false; then
sh generate_pseudotissue_sample_names.sh $downsampled_tissue_info_file $gtex_downsampled_individuals_dir $pseudotissue_sample_names_dir
fi

###################
# Above script generated this file
pseudotissue_info_file=${pseudotissue_sample_names_dir}"pseudotissue_info.txt"

###################
# Generate expression matrices for each of pseudobulk tissues
if false; then
sh generate_pseudotissue_expression_matrices.sh $pseudotissue_info_file $pseudotissue_sample_names_dir $gtex_tpm_expression_dir $gene_annotation_file $pseudotissue_expression_dir
fi


###################
# Generate covariate matrices for each of pseudobulk tissues
if false; then
sh generate_pseudotissue_covariate_matrices.sh $pseudotissue_info_file $pseudotissue_sample_names_dir $pseudotissue_expression_dir $gtex_covariate_dir $pseudotissue_covariate_dir
fi

###################
# Generate genotype matrices for each of pseudobulk tissues
if false; then
sed 1d $pseudotissue_info_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
	pseudotissue_individual_names=${pseudotissue_sample_names_dir}${pseudotissue_name}"_individual_names_plink_ready.txt"
	pseudotime_sample_names=${pseudotissue_sample_names_dir}${pseudotissue_name}"_sample_names.txt"
	sh generate_pseudotissue_genotype_matrices.sh $pseudotissue_name $pseudotissue_individual_names $pseudotime_sample_names $gtex_genotype_dir $pseudotissue_genotype_dir
done
fi




###################
# Run matrix eQTL (in each pseudotissue)
# loop through pseudotissues here
if false; then
sed 1d $pseudotissue_info_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
	sbatch matrix_eqtl_shell.sh $pseudotissue_name $pseudotissue_expression_dir $pseudotissue_genotype_dir $pseudotissue_covariate_dir $pseudotissue_eqtl_dir
done
fi


###################
# Run LMM eQTL (in each pseudotissue)
# loop through pseudotissues here (only perform for pseudotissues with sample repeat structure)
if false; then
sed 1d $pseudotissue_info_file | while read pseudotissue_name sample_size sample_repeat composit_tissue_string; do
	if [ "$sample_repeat" = "True" ]; then
		for chrom_num in $(seq 1 22); do 
			sbatch lmm_eqtl_shell.sh $pseudotissue_name $pseudotissue_expression_dir $pseudotissue_genotype_dir $pseudotissue_covariate_dir ${pseudotissue_sample_names_dir} $pseudotissue_eqtl_dir $chrom_num
		done
	fi
done
fi




###################
# Generate gene level genotype reference panels
cross_tissue_gene_list_file=${pseudotissue_expression_dir}"cross_tissue_gene_list.txt"
if false; then
for chrom_num in $(seq 1 22); do 
	sbatch generate_gene_level_genotype_reference_panel.sh $chrom_num $cross_tissue_gene_list_file $gtex_genotype_dir $genotype_reference_panel_dir
done
fi


