args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)
library(reshape2)



run_eqtl_lmm <- function(expr, geno, covariates, groups) {
	fit_full <- lmer(expr ~ geno + covariates + (1 | groups), REML=FALSE)

	# extract coefficients
	coefs <- data.frame(coef(summary(fit_full)))


	effect_size <- coefs[2,1]
	std_err <- coefs[2,2]
	t_value <- coefs[2,3]
	aggregate_pvalue <- 2 * (1 - pnorm(abs(t_value)))

	return(list(eqtl_pvalue=aggregate_pvalue, effect_size=effect_size, t_value=t_value))
}






#######################
# Command line args
#######################
matrix_eqtl_genotype_file <- args[1]
matrix_eqtl_snp_loc_file <- args[2]
matrix_eqtl_expression_file <- args[3]
matrix_eqtl_gene_loc_file <- args[4]
matrix_eqtl_covariate_file <- args[5]
matrix_eqtl_sample_repeat_file <- args[6]
cis_distance <- as.numeric(args[7])
output_file <- args[8]


##############################
# Load in data
##############################

# Load in sample repeat data
groups <- read.table(matrix_eqtl_sample_repeat_file, header=FALSE)$V1 +1

# Load in snp location and gene location data
snp_location <- read.table(matrix_eqtl_snp_loc_file, header=TRUE,sep="\t")
gene_location <- read.table(matrix_eqtl_gene_loc_file, header=TRUE, sep="\t")

# Load in covariate data
cov_df <- read.table(matrix_eqtl_covariate_file,header=TRUE,sep="\t")
cov_df <- t(as.matrix(cov_df[,2:(dim(cov_df)[2])]))

# Load in expression data
expr_df <- read.table(matrix_eqtl_expression_file,header=TRUE,sep="\t")
expr_df <- (as.matrix(expr_df[,2:(dim(expr_df)[2])]))
# To get slice: as.numeric(expr_df[row_index,])

# Load in genotype data
geno_df <- read.table(matrix_eqtl_genotype_file,header=TRUE,sep="\t")
geno_df <- (as.matrix(geno_df[,2:(dim(geno_df)[2])]))


# Open output file handle
sink(output_file)

# Print header
new_line <- paste0('SNP', "\t", 'gene' ,"\t",'beta', "\t", 't-stat', "\t", 'p-value',"\t", 'converged', "\n")
cat(new_line)




##############################
# Loop through genes
##############################
num_genes <- dim(expr_df)[1]
for (gene_index in 1:num_genes) {
	# Extract info for gene
	gene_id <- as.character(gene_location$geneid[gene_index])
	gene_chrom <- as.character(gene_location$chr[gene_index])
	gene_tss <- gene_location$s1[gene_index]
	expr_vec <- as.numeric(expr_df[gene_index,])

	# Get indices of snps in cis window of gene
	cis_window_snp_indices <- abs(snp_location$pos - gene_tss) <= cis_distance
	# subset genotype data to just cis snps
	cis_window_geno_df <- geno_df[cis_window_snp_indices,]
	cis_window_snp_loc <- snp_location[cis_window_snp_indices,]


	##############################
	# Loop through cis snps for genes
	##############################
	num_cis_snps = sum(cis_window_snp_indices)
	for (cis_snp_num in 1:num_cis_snps) {
		genotype_vec <- as.numeric(cis_window_geno_df[cis_snp_num,])
		snp_id <- as.character(cis_window_snp_loc$snp[cis_snp_num])
		snp_pos <- cis_window_snp_loc$pos[cis_snp_num]
		# Try catch (sometimes lmm fails)
		tryCatch(
		{
			lmm_results = run_eqtl_lmm(expr_vec, genotype_vec, cov_df, groups)
			new_line <- paste0(snp_id, "\t", gene_id ,"\t", lmm_results$effect_size, "\t", lmm_results$t_value, "\t", lmm_results$eqtl_pvalue, "\t", "True", "\n")
			cat(new_line)
		},
		error = function(e) {
			new_line <- paste0(snp_id, "\t", gene_id ,"\t", 0.0, "\t", 0.0, "\t", 1.0, "\t", "False", "\n")
			cat(new_line)
		}
		)
	}

}

sink()




