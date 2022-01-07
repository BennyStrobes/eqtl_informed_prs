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

input_dir <- args[1]
output_dir <- args[2]




prs_weights_file <- paste0(input_dir, "blood_white_count_prs_weights.txt")
prs_weights <- read.table(prs_weights_file, header=TRUE, sep="\t")

relative_r_squared_file <- paste0(input_dir, "blood_white_count_relative_r_squared.txt")
relative_r_squared <- read.table(relative_r_squared_file, header=TRUE, sep="\t")

corr_file <- paste0(input_dir, "blood_white_count_prs_correlation_matrix.txt")
corr_mat <- read.table(corr_file, header=TRUE, sep="\t")
corr_mat <- corr_mat[,2:(dim(corr_mat)[2])]

heatmap <- make_heatmap(corr_mat)
ggsave(heatmap, file=paste0(output_dir, "headmap.pdf"))

