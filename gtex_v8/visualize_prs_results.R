args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(reshape)

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

make_heatmap <- function(corr_mat) {
	ord <- hclust( dist(corr_mat, method = "euclidean"), method = "ward.D" )$order

	ordered_names <- colnames(corr_mat)
	row.names(corr_mat) <- ordered_names


	correlation_matrix <- as.matrix(corr_mat)
	melted_corr <- melt(correlation_matrix)

	print(melted_corr)


    melted_corr$X1 <- factor(melted_corr$X1, levels = ordered_names[ord])
    melted_corr$X2 <- factor(melted_corr$X2, levels = ordered_names[ord])
   # print(melted_corr)

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu", direction=1)
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90))
    heatmap <- heatmap + labs(x = "tissue type", y = "tissue type", fill= "Spearman Rho") + gtex_v8_figure_theme()
    return(heatmap)

}

make_prs_pc_scatter_colored_by_covariate <- function(pc1, pc2, cov, pc1_name, pc2_name, cov_name) {
    df <- data.frame(loading_1=pc1, loading_2=pc2, covariate=cov)
    plotter <- ggplot(df) + 
               geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.1) +
               gtex_v8_figure_theme() + 
               labs(x=pc1_name, y = pc2_name, color=cov_name) + 
               scale_colour_gradient2() +
               scale_color_gradient(low="pink",high="blue") +
               theme(legend.position="bottom") + 
               theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
    return(plotter)
}

make_scatter_density_plot <- function(pc1, pc2, pc1_name, pc2_name) {
    df <- data.frame(loading_1=pc1, loading_2=pc2)
    plotter <- ggplot(df, aes(x=loading_1, y=loading_2)) + 
               geom_hex(bins=70) +
                scale_fill_continuous(type = "viridis") +
                theme_bw() +
               gtex_v8_figure_theme() + 
               labs(x=pc1_name, y = pc2_name) + 
               theme(legend.position="bottom") + 
               theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
    return(plotter)

}

input_dir <- args[1]
output_dir <- args[2]




prs_weights_file <- paste0(input_dir, "blood_white_count_prs_weights.txt")
prs_weights <- read.table(prs_weights_file, header=TRUE, sep="\t")

relative_r_squared_file <- paste0(input_dir, "blood_white_count_relative_r_squared.txt")
relative_r_squared <- read.table(relative_r_squared_file, header=TRUE, sep="\t")

corr_file <- paste0(input_dir, "blood_white_count_prs_correlation_matrix.txt")
corr_mat <- read.table(corr_file, header=TRUE, sep="\t")
corr_mat <- corr_mat[,2:(dim(corr_mat)[2])]

prs_pc_file = paste0(input_dir, "blood_white_count_prs_pcs.txt")
prs_pcs = read.table(prs_pc_file, header=TRUE, sep="\t")

cov_file = paste0(input_dir, "blood_white_count_blood_covariates.txt")
cov = read.table(cov_file, header=TRUE, sep="\t")


# Make PRS PC scatter colored by covariate
output_file <- paste0(output_dir, "prs_pc1_vs_prs_pc2_scatter_colored_by_blood_white_count.pdf")
scatter <- make_prs_pc_scatter_colored_by_covariate(prs_pcs$prs_pc1, prs_pcs$prs_pc2, cov$blood_WHITE_COUNT, "PRS_PC1", "PRS_PC2", "blood_white_count")
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_pc2_vs_prs_pc3_scatter_colored_by_blood_white_count.pdf")
scatter <- make_prs_pc_scatter_colored_by_covariate(prs_pcs$prs_pc2, prs_pcs$prs_pc3, cov$blood_WHITE_COUNT, "PRS_PC2", "PRS_PC3", "blood_white_count")
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_pc1_vs_blood_white_count_scatter_colored_by_prs_pc2.pdf")
scatter <- make_prs_pc_scatter_colored_by_covariate(prs_pcs$prs_pc1, cov$blood_WHITE_COUNT, prs_pcs$prs_pc2, "PRS_PC1", "blood_white_count", "PRS_PC2")
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_pc2_vs_blood_NEUTROPHIL_PCT.pdf")
scatter <- make_scatter_density_plot(prs_pcs$prs_pc2, cov$blood_NEUTROPHIL_PCT, "PRS_PC2", "blood_neutrophil_pct")
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_pc1_vs_blood_WHITE_COUNT.pdf")
scatter <- make_scatter_density_plot(prs_pcs$prs_pc1, cov$blood_WHITE_COUNT, "PRS_PC1", "blood_white_count")
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")


#heatmap <- make_heatmap(corr_mat)
#ggsave(heatmap, file=paste0(output_dir, "headmap.pdf"))

