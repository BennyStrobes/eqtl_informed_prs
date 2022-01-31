





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

# UKBB download data
ukbb_download_data="/n/groups/price/UKBiobank/download_500K/"

# UKBB Genotype data
ukbb_genotype_data="/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/"

# UKBB Phenotype files
ukbb_pheno_file1="/n/groups/price/steven/RareVariants/Final/UKB_new_sumstats/UKB_v3.061518.tab"
ukbb_pheno_file2="/n/groups/price/UKBiobank/app19808mosaic/bloodQC/ukb4777.blood_v2.covars.tab"
ukbb_pheno_file3="/n/groups/price/UKBiobank/app10438assoc/ukb4777.processed_and_post.plinkPCs.tab.gz"


##################
# Output data
##################
# Output root directory
output_root="/n/groups/price/ben/eqtl_informed_prs/gtex_v8_tissue_specific_prs/"

# Directory containing all gtex associations (in cafeh format)
processed_gtex_associations_dir=$output_root"processed_gtex_associations/"

# Directory containing organized bivariate cafeh input data
processed_bivariate_cafeh_input_dir=$output_root"processed_bivariate_cafeh_input/"

# Directory containing bivariate CAFEH results
bivariate_cafeh_output_dir=$output_root"bivariate_cafeh_output/"

# Directory containing organized multivariate cafeh input data
processed_multivariate_cafeh_input_dir=$output_root"processed_multivariate_cafeh_input/"

# Directory containing multivariate cafeh results
multivariate_cafeh_output_dir=$output_root"multivariate_cafeh_output/"

# Directory containing filtered UKBB sample names
ukbb_sample_names_dir=$output_root"ukbb_sample_names/"

# Directory containing generated PRS
ukbb_prs_dir=$output_root"ukbb_prs/"

# Directory containing generated PRS
ukbb_multivariate_prs_dir=$output_root"ukbb_multivariate_prs/"

# Directory containing analyzed UKBB PRS results
analyzed_ukbb_prs_dir=$output_root"analyzed_ukbb_prs/"

# Directory containing analyzed UKBB PRS results
analyzed_ukbb_multivariate_prs_dir=$output_root"analyzed_ukbb_multivariate_prs/"

# Directory containing visualizations of generated
ukbb_prs_viz_dir=$output_root"visualize_ukbb_prs/"



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



trait_name="blood_WHITE_COUNT"
trait_sumstat_file=$summary_stat_dir"bolt_337K_unrelStringentBrit_MAF0.001_v3."${trait_name}".bgen.stats.gz"
sample_size="337491"
if false; then
sh process_data_for_multivariate_cafeh_analysis.sh $trait_name $trait_sumstat_file $gtex_tissue_file $genotype_reference_panel_dir $processed_gtex_associations_dir $processed_multivariate_cafeh_input_dir $sample_size $cafeh_gene_list_file
fi

version="all_tissues_.95"
if false; then
sh run_cafeh_on_multivariate_data_for_single_trait_shell.sh $trait_name $gtex_tissue_file $processed_multivariate_cafeh_input_dir $multivariate_cafeh_output_dir $version
fi










if false; then
sh generate_ukbb_sample_name_lists.sh $ukbb_download_data $ukbb_sample_names_dir
fi



if false; then
for chrom_num in $(seq 1 22); do
	sbatch generate_prs.sh $ukbb_sample_names_dir $ukbb_download_data $ukbb_genotype_data $bivariate_cafeh_output_dir $ukbb_prs_dir $chrom_num $trait_name $version
done
fi

if false; then
sh organize_prs_results.sh $ukbb_prs_dir $trait_name $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $analyzed_ukbb_prs_dir
fi

if false; then
module load R/3.5.1
model_version="beta_thresh_0.05_weighted"
Rscript visualize_prs_results.R $trait_name $bivariate_cafeh_output_dir $ukbb_prs_dir $analyzed_ukbb_prs_dir $ukbb_prs_viz_dir $model_version
fi







if false; then
for chrom_num in $(seq 1 22); do
	sbatch generate_prs.sh $ukbb_sample_names_dir $ukbb_download_data $ukbb_genotype_data $multivariate_cafeh_output_dir $ukbb_multivariate_prs_dir $chrom_num $trait_name $version
done
fi


if false; then
sh organize_prs_results.sh $ukbb_multivariate_prs_dir $trait_name $ukbb_pheno_file1 $ukbb_pheno_file2 $ukbb_pheno_file3 $analyzed_ukbb_multivariate_prs_dir
fi


if false; then
module load R/3.5.1
model_version="beta_thresh_0.05_weighted"
Rscript visualize_prs_results.R $trait_name $multivariate_cafeh_output_dir $ukbb_multivariate_prs_dir $analyzed_ukbb_multivariate_prs_dir $ukbb_prs_viz_dir $model_version
fi





