





##################
# Input data
##################
# Directory containing cafeh scripts
# NEED TO UPDATE
cafeh_code_dir="/work-zfs/abattle4/bstrober/tools/cafeh/"

# Directory containing summary statistics
# contains file called study_info explaining files and where they came from
# NEED TO UPDATE
summary_stat_dir="/n/groups/price/UKBiobank/sumstats/bolt_337K_unrelStringentBrit_MAF0.001_v3/"

#Directory that contains necessary liftover information.
##Specifically, it must contain:
#####1. 'liftOver'   --> the executable
#####2. 'hg19ToHg38.over.chain.gz'   --> for converting from hg19 to hg38
#####2. 'hg38ToHg19.over.chain.gz'   --> for converting from hg38 to hg19
liftover_directory="/n/groups/price/ben/tools/liftOver_x86/"

# File containing gtex tissues to do analysis on and their sample size
gtex_tissue_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/pseudotissue_sample_names/pseudotissue_info.txt"

# File containing list of genes used for cafeh
cafeh_gene_list_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/pseudotissue_expression/cross_tissue_gene_list.txt"

# Directory containing files of eqtl summary stats
eqtl_summary_stats_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/pseudotissue_eqtl_summary_stats/"

# Directory containing genotype reference panel data
genotype_reference_panel_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_eqtl_calling/gene_level_genotype_reference_panel/"




##################
# Output data
##################
# Output root directory
output_root="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_tissue_specific_prs/"

# Directory containing all gtex associations (in cafeh format)
processed_gtex_associations_dir=$output_root"processed_gtex_associations/"

# Directory containing all gtex associations (in cafeh format)
processed_bivariate_cafeh_input_dir=$output_root"processed_bivariate_cafeh_input/"

# Directory containing CAFEH results
bivariate_cafeh_output_dir=$output_root"bivariate_cafeh_output/"








##################
# Run analysis
##################
if false; then
sh process_gtex_associations_for_cafeh.sh $processed_gtex_associations_dir $cafeh_gene_list_file $liftover_directory $eqtl_summary_stats_dir $gtex_tissue_file $genotype_reference_panel_dir
fi

#trait_name="whr_adjusted_bmi"
#trait_sumstat_file=$summary_stat_dir"bolt_337K_unrelStringentBrit_MAF0.001_v3.body_WHRadjBMIz.bgen.stats.gz"
#sample_size="487409"

trait_name="blood_WHITE_COUNT"
trait_sumstat_file=$summary_stat_dir"bolt_337K_unrelStringentBrit_MAF0.001_v3."${trait_name}".bgen.stats.gz"
sample_size="337491"
if false; then
sh process_data_for_bivariate_cafeh_analysis.sh $trait_name $trait_sumstat_file $gtex_tissue_file $genotype_reference_panel_dir $processed_gtex_associations_dir $processed_bivariate_cafeh_input_dir $sample_size $cafeh_gene_list_file
fi

if false; then
sh run_cafeh_on_bivariate_data_for_single_trait_shell.sh $trait_name $gtex_tissue_file $processed_bivariate_cafeh_input_dir $bivariate_cafeh_output_dir
fi














# CURRENTLY HACKY, switching over to o2

##################
# Input data
##################
# UKBB download data
ukbb_download_data="/n/groups/price/UKBiobank/download_500K/"
# UKBB sumstats
ukbb_sumstats_dir="/n/groups/price/UKBiobank/sumstats/bolt_337K_unrelStringentBrit_MAF0.001_v3/"
# Beta files
cafeh_prs_betas_dir="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_tissue_specific_prs/input_data/"
# File containing estimates of number of cafeh components per tissue
num_components_per_tissue_file="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_tissue_specific_prs/input_data/cafeh_results_whr_adjusted_bmi_num_prs_components.txt"
# UKBB Phenotype files
ukbb_pheno_file1="/n/groups/price/steven/RareVariants/Final/UKB_new_sumstats/UKB_v3.061518.tab"
ukbb_pheno_file2="/n/groups/price/UKBiobank/app19808mosaic/bloodQC/ukb4777.blood_v2.covars.tab"
ukbb_pheno_file3="/n/groups/price/UKBiobank/app10438assoc/ukb4777.processed_and_post.plinkPCs.tab.gz"



##################
# Output data
##################
# Root output directory
output_root="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_tissue_specific_prs/"
# Directory containing filtered UKBB sample names
ukbb_sample_names_dir=$output_root"ukbb_sample_names/"
# Directory containing generate PRS
ukbb_prs_dir=$output_root"ukbb_prs/old/"
ukbb_prs_dir=$output_root"ukbb_prs/"
# Directory containing generate PRS
ukbb_prs_viz_dir=$output_root"visualize_ukbb_prs/"

trait_name="whr_adjusted_bmi"



if false; then
sh generate_ukbb_sample_name_lists.sh $ukbb_download_data $ukbb_sample_names_dir
fi



beta_threshold=".5"
chrom_num="1"
if false; then
sbatch generate_prs.sh $ukbb_sample_names_dir $ukbb_download_data $cafeh_prs_betas_dir $ukbb_prs_dir $chrom_num $beta_threshold $trait_name
fi

beta_thresholds=( "1" ".5" ".1" ".05" ".01" ".005")
beta_thresholds=(".05")
if false; then
for beta_threshold in "${beta_thresholds[@]}"; do	
	for chrom_num in $(seq 1 22); do
		sbatch generate_prs.sh $ukbb_sample_names_dir $ukbb_download_data $cafeh_prs_betas_dir $ukbb_prs_dir $chrom_num $beta_threshold $trait_name
	done
done
fi

if false; then
sh organize_prs_results.sh $ukbb_prs_dir $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3
fi

if false; then
module load R/3.5.1
Rscript visualize_prs_results.R $ukbb_prs_dir $ukbb_prs_viz_dir $num_components_per_tissue_file
fi








