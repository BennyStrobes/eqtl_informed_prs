args = commandArgs(trailingOnly=TRUE)
library(coloc)
library(hash)
#options(warn=1)



run_susie_coloc_for_single_gene <- function(gene_name, genotype_file, n_file, beta_file, std_err_file, tissue_names, coloc_output_dir) {
	# Load in coloc data
	beta_df <- read.table(beta_file, sep="\t", row.names=1, header=TRUE)
	std_err_df <- read.table(std_err_file, sep="\t", row.names=1, header=TRUE)
	n_df <- read.table(n_file, sep="\t", row.names=1, header=TRUE)
	geno_df <- read.table(genotype_file, header=TRUE, sep="\t", row.names=1)
	LD <- cor(geno_df)

	# Get studies observed for this gene
	observed_studies <- row.names(beta_df)
	num_studies = length(observed_studies)
	eqtl_studies <- observed_studies[1:(num_studies-1)]

	# Get snp names
	snp_names = colnames(beta_df)

	# Get Trait data in correct format for coloc
	trait_beta <- (as.numeric(beta_df[num_studies,]))
	trait_var_beta <- (as.numeric(std_err_df[num_studies,]))^2
	trait_df = data.frame(beta=trait_beta, varbeta=trait_var_beta)
	trait_data <- as.list(trait_df)
	trait_data$LD = LD
	trait_data$snp = snp_names
	trait_data$type="quant"
	trait_data$sdY = 1.0
	trait_data$N = n_df[num_studies,1]

	# Run susie on trait data
	trait_data_susie = runsusie(trait_data)
	num_trait_components = length(trait_data_susie[16]$sets$cs)

	if (num_trait_components > 0) {
		# Now loop through tissues
		for (eqtl_study_num in 1:length(eqtl_studies)) {
			eqtl_study_name = eqtl_studies[eqtl_study_num]
			print(eqtl_study_name)
			# Get Trait data in correct format for coloc
			eqtl_beta <- (as.numeric(beta_df[eqtl_study_num,]))
			eqtl_var_beta <- (as.numeric(std_err_df[eqtl_study_num,]))^2
			eqtl_df = data.frame(beta=eqtl_beta, varbeta=eqtl_var_beta)
			eqtl_data <- as.list(eqtl_df)
			eqtl_data$LD = LD
			eqtl_data$snp = snp_names
			eqtl_data$type="quant"
			eqtl_data$sdY = 1.0
			eqtl_data$N = n_df[eqtl_study_num,1]	

			# Run susie on trait data
			eqtl_data_susie = runsusie(eqtl_data)
			num_eqtl_components = length(eqtl_data_susie[16]$sets$cs)
			if (num_eqtl_components > 0) {
				# Test for coloc
				res=coloc.susie(trait_data_susie,eqtl_data_susie)
				if (max(res$summary[,8]) > .5) {
					print("HIT")

					saveRDS(res, file=paste0(gene_name,"_", eqtl_study_name, "_coloc.RDS"))
					saveRDS(eqtl_data_susie, file=paste0(gene_name,"_", eqtl_study_name, "_eqtl_susie.RDS"))
					saveRDS(trait_data_susie, file=paste0(gene_name,"_", eqtl_study_name, "_trait_susie.RDS"))
				}
			}
		}


	}



}


run_coloc_for_single_gene <- function(gene_name, genotype_file, n_file, beta_file, std_err_file, tissue_names, tissue_name_to_position, coloc_output_dir, trait_name) {
	# Coloc thresholds
	coloc_thresholds <- c(.5, .7, .9, .95, .99)

	# Load in coloc data
	beta_df <- read.table(beta_file, sep="\t", row.names=1, header=TRUE)
	std_err_df <- read.table(std_err_file, sep="\t", row.names=1, header=TRUE)
	n_df <- read.table(n_file, sep="\t", row.names=1, header=TRUE)
	#geno_df <- read.table(genotype_file, header=TRUE, sep="\t", row.names=1)
	#LD <- cor(geno_df)

	# Get studies observed for this gene
	observed_studies <- row.names(beta_df)
	num_studies = length(observed_studies)
	eqtl_studies <- observed_studies[1:(num_studies-1)]

	# Get snp names
	snp_names = colnames(beta_df)

	# Get Trait data in correct format for coloc
	trait_beta <- (as.numeric(beta_df[num_studies,]))
	trait_var_beta <- (as.numeric(std_err_df[num_studies,]))^2
	trait_df = data.frame(beta=trait_beta, varbeta=trait_var_beta)
	trait_data <- as.list(trait_df)
	#trait_data$LD = LD
	trait_data$snp = snp_names
	trait_data$type="quant"
	trait_data$sdY = 1.0
	trait_data$N = n_df[num_studies,1]

	# Initilize output
	predicted_effects_list = list()
	coloc_at_threshold_arr <- c()
	for (threshold_iter in 1:length(coloc_thresholds)) {
		predicted_effects_list[[threshold_iter]] = matrix(0, length(snp_names), length(tissue_names))
		coloc_at_threshold_arr <- c(coloc_at_threshold_arr, FALSE)
	}


	pp_across_tissues = matrix(0, length(tissue_names), 5)
	pp_across_tissues[,1] = 1.0

	# Now loop through tissues
	for (eqtl_study_num in 1:length(eqtl_studies)) {
		eqtl_study_name = eqtl_studies[eqtl_study_num]
		tissue_position = tissue_name_to_position[[eqtl_study_name]]
		# Get Trait data in correct format for coloc
		eqtl_beta <- (as.numeric(beta_df[eqtl_study_num,]))
		eqtl_var_beta <- (as.numeric(std_err_df[eqtl_study_num,]))^2
		eqtl_df = data.frame(beta=eqtl_beta, varbeta=eqtl_var_beta)
		eqtl_data <- as.list(eqtl_df)
		eqtl_data$snp = snp_names
		eqtl_data$type="quant"
		eqtl_data$sdY = 1.0
		eqtl_data$N = n_df[eqtl_study_num,1]	

		# Test for coloc
		res <- coloc.abf(dataset1=eqtl_data, dataset2=trait_data)

		pph4 <- res$summary[6]

		pp_across_tissues[tissue_position,] = res$summary[2:length(res$summary)]

		snp_probs = res$results$SNP.PP.H4
		snp_probs[snp_probs<.01] = 0

		for (threshold_iter in 1:length(coloc_thresholds)) {
			coloc_threshold = coloc_thresholds[threshold_iter]
			if (pph4 > coloc_threshold) {
				coloc_at_threshold_arr[threshold_iter] = TRUE
				predicted_effects_list[[threshold_iter]][,tissue_position] = snp_probs*trait_beta*pph4
			}
		}
	}
	pp_across_tissues_file <- paste0(coloc_output_dir, trait_name, "_", gene_name, "_coloc_posterior_probabilities.txt")
	write.table(pp_across_tissues, file=pp_across_tissues_file, row.names=tissue_names, col.names=c("PPH0", "PPH1", "PPH2", "PPH3", "PPH4"), sep="\t", quote = FALSE)


	for (threshold_iter in 1:length(coloc_thresholds)) {
		coloc_threshold = coloc_thresholds[threshold_iter]
		if (coloc_at_threshold_arr[threshold_iter]) {
			predicted_effect_size_file <- paste0(coloc_output_dir, trait_name, "_", gene_name, "_coloc_", coloc_threshold, "_predicted_effect_sizes.txt")
			write.table(predicted_effects_list[[threshold_iter]], file=predicted_effect_size_file, row.names=snp_names, col.names=tissue_names, sep="\t", quote = FALSE)
		}
	}



}



#####################
# Command line args
#####################
gene_file = args[1]
trait_name = args[2]
gtex_tissue_file = args[3]
coloc_output_dir = args[4]
job_number = as.numeric(args[5])
total_jobs = as.numeric(args[6])

# Load in tissue names
tissue_df <- read.table(gtex_tissue_file, header=TRUE, sep="\t")
tissue_names <- tissue_df$pseudotissue_name
tissue_name_to_position <- hash()
for (tissue_iter in 1:length(tissue_names)) {
	tissue_name_to_position[[tissue_names[tissue_iter]]] = tissue_iter
}


# Load in gene data frame
gene_df <- read.table(gene_file,header=TRUE, sep="\t")
# Get total number of genes
num_genes <- dim(gene_df)[1]


# For parallelization purposes, determine which genes to test in this thread
tasks_per_job = floor(num_genes/total_jobs) + 1
start_task = floor(job_number*tasks_per_job + 1)
end_task = floor((job_number + 1)*tasks_per_job)
if (end_task > num_genes) {
	end_task = num_genes
}



# Loop through all genes in this thread
for (gene_num in start_task:end_task) {
	print(gene_num)
	gene_name <- gene_df$gene_id[gene_num]
	genotype_file <- gene_df$genotype_file[gene_num]
	n_file <- gene_df$n_file[gene_num]
	beta_file <- gene_df$beta_file[gene_num]
	std_err_file <- gene_df$std_err_file[gene_num]
	tryCatch(
		{
			#run_susie_coloc_for_single_gene(gene_name, genotype_file, n_file, beta_file, std_err_file, tissue_names, coloc_output_dir)
			run_coloc_for_single_gene(gene_name, genotype_file, n_file, beta_file, std_err_file, tissue_names, tissue_name_to_position, coloc_output_dir, trait_name)

		},
			error = function(e) {
				print('ERRORR')
		}
	)

}
