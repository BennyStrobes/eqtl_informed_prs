args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(reshape)

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}
gtex_v8_figure_theme_small <- function() {
	return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),axis.text=element_text(size=8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=8), legend.title = element_text(size=8)))
}

make_number_of_coloc_components_for_single_method_seperate_plots <- function(coloc_output_dir, method, trait_name, coloc_thresholds) {
	# First make global df keeping track of number of colocalization across thresholds
	plot_arr <- list()
	for (coloc_thresh_iter in 1:length(coloc_thresholds)) {
		num_coloc_arr <- c()
		tissue_arr <- c()
		coloc_threshold_arr <- c()

		coloc_threshold <- coloc_thresholds[coloc_thresh_iter]
		num_components_per_tissue_file <- paste0(coloc_output_dir, method, "_results_", trait_name, "_", coloc_threshold, "_num_prs_components.txt")
		tissue_df = read.table(num_components_per_tissue_file, header=TRUE, sep="\t")

		num_coloc_arr <- c(num_coloc_arr, tissue_df$num_cafeh_components)
		tissue_arr <- c(tissue_arr, as.character(tissue_df$tissue))
		coloc_threshold_arr <- c(coloc_threshold_arr, rep(coloc_threshold, dim(tissue_df)[1]))
		tissue_names <- as.character(tissue_df$tissue)

		df <- data.frame(num_coloc=num_coloc_arr, tissue=factor(tissue_arr, levels=tissue_names))

  		p <- ggplot(data=df, aes(x=tissue, y=num_coloc)) +
    		geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    		gtex_v8_figure_theme() +
     		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    		theme(axis.text.y = element_text(size=8)) +
    		theme(axis.title = element_text(size=8)) +
    		labs(y=paste0("Num coloc \n(", coloc_threshold,")"), x="")


     	if (coloc_thresh_iter != length(coloc_thresholds)) {
      		plot_arr[[coloc_thresh_iter]] <- p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  		} else {
    		plot_arr[[coloc_thresh_iter]] <- p
  		}

	}

	merged = plot_grid(plotlist=plot_arr, ncol=1, rel_heights=c(1,1, 1, 1,2.8))
	return(merged)


}


make_mixture_model_causal_prior_probability <- function(coloc_output_dir, trait_name, version) {
	file_name <- paste0(coloc_output_dir, trait_name, "_coloc_mm_", version, "_causal_prior_probability.txt")
	df <- read.table(file_name, header=TRUE, sep="\t")
  p <- ggplot(data=df, aes(x=tissue_name, y=prior_probability)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
     theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) +
    labs(y="Mixture model weights", x="") 
   return(p)

}
make_number_of_coloc_components_for_single_method <- function(coloc_output_dir, method, trait_name, coloc_thresholds) {
	# First make global df keeping track of number of colocalization across thresholds
	num_coloc_arr <- c()
	tissue_arr <- c()
	coloc_threshold_arr <- c()

	for (coloc_thresh_iter in 1:length(coloc_thresholds)) {
		coloc_threshold <- coloc_thresholds[coloc_thresh_iter]
		num_components_per_tissue_file <- paste0(coloc_output_dir, method, "_results_", trait_name, "_", coloc_threshold, "_num_prs_components.txt")
		tissue_df = read.table(num_components_per_tissue_file, header=TRUE, sep="\t")

		num_coloc_arr <- c(num_coloc_arr, tissue_df$num_cafeh_components)
		tissue_arr <- c(tissue_arr, as.character(tissue_df$tissue))
		coloc_threshold_arr <- c(coloc_threshold_arr, rep(coloc_threshold, dim(tissue_df)[1]))
		tissue_names <- as.character(tissue_df$tissue)
	}

	df <- data.frame(num_coloc=num_coloc_arr, tissue=factor(tissue_arr, levels=tissue_names), coloc_threshold=coloc_threshold_arr)

  p <- ggplot(data=df, aes(x=tissue, y=num_coloc, fill=coloc_threshold)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
     theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) +
    labs(y="Number of colocalizations", x="", fill="Coloc\nthreshold") 
  return(p) 

}


make_histogram_of_pph4s <- function(pph4_vec, tissue_name, pph4_threshold) {

	filtered_vec <- pph4_vec[pph4_vec > pph4_threshold]


	df = data.frame(pph4=filtered_vec)

	p <- ggplot(df, aes(x=pph4)) +
 	 geom_histogram(color="darkblue", fill="lightblue") + 
 	 gtex_v8_figure_theme() +
 	 labs(x=paste0(tissue_name, " pph4"))
 	return(p)

}

make_histogram_of_overlapping_pph4s <- function(pph4_df, pph4_threshold) {
	pph4s <- c()
	tissue_names <- c()

	tissue_pph4 <- pph4_df$Whole_Blood[pph4_df$Whole_Blood > 0.1]
	pph4s <- c(pph4s, tissue_pph4)
	tissue_names <- c(tissue_names, rep("Whole_Blood", length(tissue_pph4)))

	tissue_pph4 <- pph4_df$Cells_Cultured_fibroblasts[pph4_df$Cells_Cultured_fibroblasts > 0.1]
	pph4s <- c(pph4s, tissue_pph4)
	tissue_names <- c(tissue_names, rep("Fibroblasts", length(tissue_pph4)))

	df <- data.frame(pph4=pph4s, tissue=factor(tissue_names))

	p <- ggplot(df, aes(x=pph4, color=tissue)) +
  		geom_histogram(fill="white", alpha=0.5, position="identity") +
  		gtex_v8_figure_theme() + 
  		theme(legend.position="bottom") + 
  		labs(x="pph4")
  	return(p)

}

causal_pph4_scatter <- function(coloc_output_dir, trait_name, tissue_name) {
	pphc_matrix_v1_file <- paste0(coloc_output_dir, trait_name, "_pphc_v1_matrix.txt")
	pphc_matrix_v2_file <- paste0(coloc_output_dir, trait_name, "_pphc_v2_matrix.txt")
	pph4_matrix_file <- paste0(coloc_output_dir, trait_name, "_pph4_matrix.txt")

	pphc_matrix_v1 <- read.table(pphc_matrix_v1_file, header=TRUE,sep="\t")
	pphc_matrix_v2 <- read.table(pphc_matrix_v2_file, header=TRUE,sep="\t")
	pph4_matrix <- read.table(pph4_matrix_file, header=TRUE,sep="\t")


	df <- data.frame(p_causal_v1=pphc_matrix_v1[[tissue_name]], p_causal_v2=pphc_matrix_v2[[tissue_name]], pph4=pph4_matrix[[tissue_name]])

	df = df[df$pph4 > .5,]
	p <- ggplot(df, aes(x=p_causal_v1, y=p_causal_v2, colour=pph4)) + geom_point() +
	gtex_v8_figure_theme()

	return(p)
}

distribution_of_causal_probabilities_histogram <- function(coloc_output_dir, trait_name, causal_suffix, pph4_threshold) {
	pphc_matrix_file <- paste0(coloc_output_dir, trait_name, causal_suffix)
	pph4_matrix_file <- paste0(coloc_output_dir, trait_name, "_pph4_matrix.txt")
	pph4_matrix <- read.table(pph4_matrix_file, header=TRUE,sep="\t")
	pphc_matrix <- read.table(pphc_matrix_file, header=TRUE,sep="\t")
	tissue_names <- as.character(colnames(pph4_matrix)[2:length(colnames(pph4_matrix))])


	plot_arr <- list()
	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		indices <- pph4_matrix[[tissue_name]] > pph4_threshold

		df <- data.frame(causal_prob_arr=pphc_matrix[[tissue_name]][indices])
		print(tissue_iter)
		p <- ggplot(df, aes(x=causal_prob_arr)) +
 	 		geom_histogram(color="darkblue", fill="lightblue") + 
 	 		gtex_v8_figure_theme_small() +
 	 		labs(x=paste0(tissue_name, " causal probability"))
 	 	plot_arr[[tissue_iter]] = p

	}
	joint <- plot_grid(plotlist=plot_arr, ncol=2)

  return(joint) 	

}

distribution_of_components_causal_bar_plot <- function(coloc_output_dir, trait_name, causal_suffix) {
	pphc_matrix_file <- paste0(coloc_output_dir, trait_name, causal_suffix)
	pph4_matrix_file <- paste0(coloc_output_dir, trait_name, "_pph4_matrix.txt")
	pph4_matrix <- read.table(pph4_matrix_file, header=TRUE,sep="\t")
	pphc_matrix <- read.table(pphc_matrix_file, header=TRUE,sep="\t")
	tissue_names <- as.character(colnames(pph4_matrix)[2:length(colnames(pph4_matrix))])

	thresholds <- c(.5, .7, .9)
	causal_prob_arr <- c()
	tissue_arr <- c()
	threshold_arr <- c()
	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		for (threshold_iter in 1:length(thresholds)) {
			pph4_threshold <- thresholds[threshold_iter]
			indices <- pph4_matrix[[tissue_name]] > pph4_threshold
			#fraction <- sum(pphc_matrix[[tissue_name]][indices])/sum(indices)

			causal_prob_arr <- c(causal_prob_arr, pphc_matrix[[tissue_name]][indices])
			N <- length(pphc_matrix[[tissue_name]][indices])
			tissue_arr <- c(tissue_arr, rep(tissue_name, N))
			threshold_arr <- c(threshold_arr, rep(pph4_threshold, N))
		}

	}
	df <- data.frame(tissue=factor(tissue_arr, levels=tissue_names), causal_prob_arr=causal_prob_arr, pph4_threshold=factor(threshold_arr, levels=thresholds))
	  p <- ggplot(data=df, aes(x=tissue, y=causal_prob_arr, fill=pph4_threshold)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
     theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) +
    labs(y="Causal probability", x="", fill="pph4\nthreshold") 
  return(p) 	

}

fraction_of_components_causal_bar_plot <- function(coloc_output_dir, trait_name, causal_suffix) {
	pphc_matrix_file <- paste0(coloc_output_dir, trait_name, causal_suffix)
	pph4_matrix_file <- paste0(coloc_output_dir, trait_name, "_pph4_matrix.txt")
	pph4_matrix <- read.table(pph4_matrix_file, header=TRUE,sep="\t")
	pphc_matrix <- read.table(pphc_matrix_file, header=TRUE,sep="\t")
	tissue_names <- as.character(colnames(pph4_matrix)[2:length(colnames(pph4_matrix))])

	thresholds <- c(.1, .3, .5, .7, .9)
	fraction_arr <- c()
	tissue_arr <- c()
	threshold_arr <- c()
	for (tissue_iter in 1:length(tissue_names)) {
		tissue_name <- tissue_names[tissue_iter]
		for (threshold_iter in 1:length(thresholds)) {
			pph4_threshold <- thresholds[threshold_iter]
			indices <- pph4_matrix[[tissue_name]] > pph4_threshold
			fraction <- sum(pphc_matrix[[tissue_name]][indices])/sum(indices)

			fraction_arr <- c(fraction_arr, fraction)
			tissue_arr <- c(tissue_arr, tissue_name)
			threshold_arr <- c(threshold_arr, pph4_threshold)
		}

	}
	df <- data.frame(tissue=factor(tissue_arr, levels=tissue_names), causal_fraction=fraction_arr, pph4_threshold=factor(threshold_arr, levels=thresholds))
	  p <- ggplot(data=df, aes(x=tissue, y=causal_fraction, fill=pph4_threshold)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
     theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) +
    labs(y="fraction of colocalizations causal", x="", fill="pph4\nthreshold") 
  return(p) 	

}

effective_number_of_components_bar_plot <- function(coloc_output_dir, trait_name, causal_suffix) {
	pphc_matrix_file <- paste0(coloc_output_dir, trait_name, causal_suffix)
	pph4_matrix_file <- paste0(coloc_output_dir, trait_name, "_pph4_matrix.txt")
	pph4_matrix <- read.table(pph4_matrix_file, header=TRUE,sep="\t")
	pphc_matrix <- read.table(pphc_matrix_file, header=TRUE,sep="\t")
	pph4_thresholds <- c(0.0, .1, .3, .5, .7, .9)

	tissue_names <- as.character(colnames(pph4_matrix)[2:length(colnames(pph4_matrix))])


	pph4_max <- apply(pph4_matrix[, 2:(dim(pph4_matrix)[2])], 1, max)

	tissue_arr <- c()
	threshold_arr <- c()
	effective_num_arr <- c()

	for (pph4_threshold_iter in 1:length(pph4_thresholds)) {
		pph4_threshold <- pph4_thresholds[pph4_threshold_iter]

		indices <- pph4_max > pph4_threshold

		for (tissue_iter in 1:length(tissue_names)) {
			tissue_name <- tissue_names[tissue_iter]

			effective_num <- sum((pphc_matrix[[tissue_name]][indices])*(pph4_matrix[[tissue_name]][indices]))
			effective_num_arr <- c(effective_num_arr, effective_num)
			tissue_arr <- c(tissue_arr, tissue_name)
			threshold_arr <- c(threshold_arr, pph4_threshold)
		}

	}
	df <- data.frame(tissue=factor(tissue_arr, levels=tissue_names), threshold=factor(threshold_arr, levels=pph4_thresholds), effective_number=effective_num_arr)

	  p <- ggplot(data=df, aes(x=tissue, y=effective_number, fill=threshold)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
     theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) +
    labs(y="Effective number of\n causal colocalizations", x="", fill="pph4\nthreshold") 
  return(p) 
}



coloc_output_dir <- args[1]
visualize_coloc_dir <- args[2]
trait_name <- args[3]


output_dir <- paste0(visualize_coloc_dir, trait_name, "_")



if (FALSE) {

###################################
# Histogram of PPH4s in each tissue (side by side)
###################################
pph4_mat_file <- paste0(coloc_output_dir, trait_name, "_pph4_matrix.txt")
pph4_mat <- read.table(pph4_mat_file, header=TRUE, sep="\t")

histo_wb <- make_histogram_of_pph4s(pph4_mat$Whole_Blood, "Whole_Blood", 0.5)
histo_fibro <- make_histogram_of_pph4s(pph4_mat$Cells_Cultured_fibroblasts, "Fibroblasts", 0.5)
joint_histo <- plot_grid(histo_wb, histo_fibro, ncol=1)

output_file <- paste0(output_dir, "_pph4_histogram.pdf")
ggsave(joint_histo, file=output_file, width=7.2, height=5.0, units="in")


###################################
# Histogram of PPH4s in each tissue (on top of each other)
###################################
joint_histo <- make_histogram_of_overlapping_pph4s(pph4_mat, 0.5)
output_file <- paste0(output_dir, "_pph4_overlapping_histogram.pdf")
ggsave(joint_histo, file=output_file, width=7.2, height=5.0, units="in")
}

###################################
# Calculate distribution of components that are causal in each tissue
###################################
pph4_threshold <- .5
histo <- distribution_of_causal_probabilities_histogram(coloc_output_dir, trait_name, "_pphc_v1_matrix.txt", pph4_threshold)
output_file <- paste0(output_dir, "_causal_v1_distribution_of_", pph4_threshold, "_causal_probabilities_histogram.pdf")
ggsave(histo, file=output_file, width=7.2, height=10.0, units="in")

pph4_threshold <- .9
histo <- distribution_of_causal_probabilities_histogram(coloc_output_dir, trait_name, "_pphc_v1_matrix.txt", pph4_threshold)
output_file <- paste0(output_dir, "_causal_v1_distribution_of_", pph4_threshold, "_causal_probabilities_histogram.pdf")
ggsave(histo, file=output_file, width=7.2, height=10.0, units="in")


###################################
# Compare two versions of causal PPH4s with scatter plot
###################################
tissue_name <- "Whole_Blood"
scatter <- causal_pph4_scatter(coloc_output_dir, trait_name, tissue_name)
output_file <- paste0(output_dir, "_causal_pph4_scatter_in_", tissue_name, ".pdf")
ggsave(scatter, file=output_file, width=7.2, height=5.0, units="in")

tissue_name <- "Cells_Cultured_fibroblasts"
scatter <- causal_pph4_scatter(coloc_output_dir, trait_name, tissue_name)
output_file <- paste0(output_dir, "_causal_pph4_scatter_in_", tissue_name, ".pdf")
ggsave(scatter, file=output_file, width=7.2, height=5.0, units="in")


###################################
# Calculate effective number of causal components bar plot
###################################
barplot <- effective_number_of_components_bar_plot(coloc_output_dir, trait_name, "_pphc_v1_matrix.txt")
output_file <- paste0(output_dir, "_causal_v1_effective_number_of_components_barplot.pdf")
ggsave(barplot, file=output_file, width=7.2, height=5.0, units="in")

barplot <- effective_number_of_components_bar_plot(coloc_output_dir, trait_name, "_pphc_v2_matrix.txt")
output_file <- paste0(output_dir, "_causal_v2_effective_number_of_components_barplot.pdf")
ggsave(barplot, file=output_file, width=7.2, height=5.0, units="in")


###################################
# Calculate fraction of components that are causal in each tissue
###################################
barplot <- fraction_of_components_causal_bar_plot(coloc_output_dir, trait_name, "_pphc_v1_matrix.txt")
output_file <- paste0(output_dir, "_causal_v1_fraction_of_components_causal_barplot.pdf")
ggsave(barplot, file=output_file, width=7.2, height=5.0, units="in")


barplot <- fraction_of_components_causal_bar_plot(coloc_output_dir, trait_name, "_pphc_v2_matrix.txt")
output_file <- paste0(output_dir, "_causal_v2_fraction_of_components_causal_barplot.pdf")
ggsave(barplot, file=output_file, width=7.2, height=5.0, units="in")


###################################
# Calculate distribution of components that are causal in each tissue
###################################
barplot <- distribution_of_components_causal_bar_plot(coloc_output_dir, trait_name, "_pphc_v1_matrix.txt")
output_file <- paste0(output_dir, "_causal_v1_distribution_of_components_causal_boxplot.pdf")
ggsave(barplot, file=output_file, width=7.2, height=5.0, units="in")


barplot <- distribution_of_components_causal_bar_plot(coloc_output_dir, trait_name, "_pphc_v2_matrix.txt")
output_file <- paste0(output_dir, "_causal_v2_distribution_of_components_causal_boxplot.pdf")
ggsave(barplot, file=output_file, width=7.2, height=5.0, units="in")



##################

methods = c("coloc", "causal_v1_coloc", "causal_v2_coloc", "mm_v1_coloc", "mm_v2_coloc") 
coloc_thresholds <- c("0.1", "0.3", "0.5", "0.7", "0.9")

###################################
# mixture model causal prior probability
###################################
barplot <- make_mixture_model_causal_prior_probability(coloc_output_dir, trait_name, "v1")
output_file <- paste0(output_dir, "_mixture_model_v1_causal_prior_prob.pdf")
ggsave(barplot, file=output_file, width=7.2, height=4.0, units="in")

barplot <- make_mixture_model_causal_prior_probability(coloc_output_dir, trait_name, "v2")
output_file <- paste0(output_dir, "_mixture_model_v2_causal_prior_prob.pdf")
ggsave(barplot, file=output_file, width=7.2, height=4.0, units="in")

###################################
# Num coloc components per tissue bar plot
# Seperate color for each colocalization threshold
###################################
for (method_iter in 1:length(methods)) {
	method = methods[method_iter]

	barplot <- make_number_of_coloc_components_for_single_method(coloc_output_dir, method, trait_name, coloc_thresholds)
	output_file <- paste0(output_dir, method, "_number_of_colocalization_bar_plot.pdf")
	ggsave(barplot, file=output_file, width=7.2, height=4.0, units="in")

}


###################################
# Num coloc components per tissue bar plot
# Seperate plot for each colocalization threshold
###################################

for (method_iter in 1:length(methods)) {
	method = methods[method_iter]
	barplot <- make_number_of_coloc_components_for_single_method_seperate_plots(coloc_output_dir, method, trait_name, coloc_thresholds)
	output_file <- paste0(output_dir, method, "_number_of_colocalization_seperate_bar_plot.pdf")
	ggsave(barplot, file=output_file, width=7.2, height=8.0, units="in")

}
