args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(reshape)

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

make_prs_pca_principal_component_heatmap <- function(matty, ordered_tissues) {
  tissue_names <- as.character(matty$sample_name)
  matty <- matty[,2:(dim(matty)[2])]
  row.names(matty) <- tissue_names
  ordered_components = as.character(colnames(matty))



  correlation_matrix <- as.matrix(matty)
  melted_corr <- melt(correlation_matrix)
  melted_corr$X1 <- factor(melted_corr$X1, levels = ordered_tissues)
  melted_corr$X2 <- factor(melted_corr$X2, levels = ordered_components)

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    heatmap <- heatmap + scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000") 
    heatmap <- heatmap + theme(text = element_text(size=18), panel.background = element_blank(), axis.text.x = element_text(angle = 90))
    heatmap <- heatmap + labs(x = "Tissues", y = "PRS Principal Components", fill= "") + gtex_v8_figure_theme()
    return(heatmap)


}


make_heatmap <- function(corr_mat) {

  num_dimensions <- dim(corr_mat)[1]
  corrz <- c()
  for (dimension_num in 1:num_dimensions) {
    for (dimension_num2 in (dimension_num):num_dimensions) {
      if (dimension_num != dimension_num2) {
        corrz <- c(corrz, corr_mat[dimension_num, dimension_num2])
      }
    }

  }
 maxy = max(corrz)

 for (dimension_num in 1:num_dimensions) {
  corr_mat[dimension_num, dimension_num] = maxy
 }



	ord <- hclust( dist(corr_mat, method = "euclidean"), method = "ward.D" )$order

	ordered_names <- colnames(corr_mat)
	row.names(corr_mat) <- ordered_names


	correlation_matrix <- as.matrix(corr_mat)
	melted_corr <- melt(correlation_matrix)



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

make_cafeh_eqtl_component_vs_cafeh_coloc_component_scatter <- function(tissue_df) {
    plotter <- ggplot(tissue_df) + 
               geom_point(aes(x=num_cafeh_eqtl_components, y=num_cafeh_components)) +
               gtex_v8_figure_theme() + 
               labs(x="Number of eQTL CAFEH components", y = "Number of colocalizing CAFEH components") + 
               theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
    return(plotter)
}

make_tissue_eqtl_component_expressed_genes_scatter <- function(tissue_df) {
    plotter <- ggplot(tissue_df) + 
               geom_point(aes(x=num_eqtl_components, y=avg_num_expressed_genes)) +
               gtex_v8_figure_theme() + 
               labs(x="CAFEH eQTL components", y = "Num expressed genes") + 
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


make_pairwise_prs_correlation_histogram <- function(corr_mat) {
  num_dimensions <- dim(corr_mat)[1]
  corrz <- c()
  for (dimension_num in 1:num_dimensions) {
    for (dimension_num2 in (dimension_num):num_dimensions) {
      if (dimension_num != dimension_num2) {
        corrz <- c(corrz, corr_mat[dimension_num, dimension_num2])
      }
    }

  }
  df <- data.frame(corrz=corrz)

  p <- ggplot(df, aes(x=corrz))+
    geom_histogram(color="darkblue", fill="lightblue") +
    gtex_v8_figure_theme() + 
    labs(x="Pairwise tissue spearman correlation") 
  return(p)

}

make_r_squared_bar_plot_with_standard_errors <- function(relative_r_squared, ordered_studies) {
  relative_r_squared$prs_name = factor(relative_r_squared$prs_name, levels=ordered_studies)
    p <- ggplot(data=relative_r_squared, aes(x=prs_name, y=relative_r_squared)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(legend.position="top") +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    labs(y="PRS Relative R^2", x="", fill="") +
    geom_errorbar(aes(ymin=relative_r_squared-(relative_r_squared_standard_error), ymax=relative_r_squared+(relative_r_squared_standard_error)), position = position_dodge(), width = .75, size=.2)

  return(p)
}

make_prs_weight_bar_plot_with_standard_errors <- function(prs_weights, p_active_threshold, y_axis_label) {
  prs_weights$prs_name = factor(prs_weights$prs_name)
    p <- ggplot(data=prs_weights, aes(x=prs_name, y=weight)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(legend.position="top") +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) +
    labs(y=y_axis_label, x="", fill="", title=p_active_threshold) +
    geom_errorbar(aes(ymin=weight-(weight_standard_error), ymax=weight+(weight_standard_error)), position = position_dodge(), width = .75, size=.2)

  return(p)

}

make_tcsc_h2_barplot_with_standard_errors <- function(df) {
    p <- ggplot(data=df, aes(x=Tissue, y=tcsc_h2)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(legend.position="top") +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) +
    labs(y="TCSC h2", x="", fill="") +
    geom_errorbar(aes(ymin=tcsc_h2-(tcsc_h2_se), ymax=tcsc_h2+(tcsc_h2_se)), position = position_dodge(), width = .75, size=.2)

  return(p)
}

make_prs_weight_per_eqtl_component_bar_plot_with_standard_errors <- function(prs_weights, ordered_studies) {
  prs_weights$prs_name = factor(prs_weights$prs_name, levels=ordered_studies)
    p <- ggplot(data=prs_weights, aes(x=prs_name, y=weight)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(legend.position="top") +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    labs(y="PRS NNLS weight/eqtl component", x="", fill="") +
    geom_errorbar(aes(ymin=weight-(weight_standard_error), ymax=weight+(weight_standard_error)), position = position_dodge(), width = .75, size=.2)

  return(p)

}

make_fraction_of_eqtl_components_that_colocalize <- function(tissue_df) {
  tissue_df$fraction = tissue_df$num_cafeh_components/tissue_df$num_cafeh_eqtl_components
  print(head(tissue_df))
  p <- ggplot(data=tissue_df, aes(x=tissue, y=fraction)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    labs(y="Fraction of eQTL components that coloc", x="", fill="")
  return(p)
}

make_num_coloc_component_per_eqtl_component_bar_plot <- function(tissue_df, title) {
  tissue_df$fraction = tissue_df$num_cafeh_components/tissue_df$num_eqtl_components
  p <- ggplot(data=tissue_df, aes(x=tissue, y=fraction)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    labs(y="Number of colocalizing components/eqtl components", x="", fill="",title=title) +
        theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) 
  return(p) 
}

make_num_coloc_component_bar_plot <- function(tissue_df, title) {
  #tissue_df$fraction = tissue_df$num_cafeh_components/tissue_df$num_cafeh_eqtl_components
  p <- ggplot(data=tissue_df, aes(x=tissue, y=num_cafeh_components)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    labs(y="Number of colocalizing components", x="", fill="", title=title) +
        theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) 
  return(p) 
}

make_learned_coloc_prior_barplot <- function(df, trait_name) {
   p <- ggplot(data=df, aes(x=V1, y=V5)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    labs(y="P4", x="", fill="", title=trait_name) +
        theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=8)) +
    theme(axis.text.y = element_text(size=8)) +
    theme(axis.title = element_text(size=8)) 
  return(p)  
}

prs_principal_component_number_of_components_scatter <- function(pc_vec, tissue_names_pc, tissue_df, component_num) {
  pcs <- c()
  number_of_components <- c()
  ordered_tissues <- c()

  for (tissue_iter in 1:length(tissue_names_pc)) {
    tissue_name <- as.character(tissue_names_pc[tissue_iter])
    tissue_pc <- pc_vec[tissue_iter]
    tissue_num_components <- tissue_df[tissue_df$tissue==tissue_name,]$num_cafeh_components[1]

    pcs <- c(pcs, tissue_pc)
    number_of_components <- c(number_of_components, tissue_num_components)
    ordered_tissues <- c(ordered_tissues, tissue_name)
  }

  df <- data.frame(number_of_components=number_of_components, prs_pc=pcs, tissue=factor(ordered_tissues))
     plotter <- ggplot(df) + 
               geom_point(aes(x=number_of_components, y=prs_pc)) +
               gtex_v8_figure_theme() + 
               labs(x="Number of CAFEH components", y = paste0("PRS Principal Component ", component_num)) + 
               theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) 

    return(plotter)

}

make_pc_variance_explained_line_plot <- function(variance_explained) {
  num_pcs = length(variance_explained$variance_explained)
  df <- data.frame(ve = variance_explained$variance_explained, pc_num=1:num_pcs)

  # PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=ve)) +
                geom_line() +
                geom_point() +
                ylim(0,max(variance_explained$variance_explained) + .002) + 
                scale_x_continuous(breaks=seq(0,(num_pcs),1)) +
                labs(x = "PRS PC", y = "PVE") + 
                gtex_v8_figure_theme() 

    return(line_plot)
}


temp_correlation <- function(vec, expr_pcs) {
  expr_pcs <- expr_pcs[,2:(dim(expr_pcs)[2])]

  # Get stdev of each column of expr_pcs and remove columns with no variation
  valid_columns = c()
  for (col_num in 1:(dim(expr_pcs)[2])) {
    vary = var(expr_pcs[,col_num], na.rm=TRUE)
    if (vary != 0.0) {
      valid_columns <- c(valid_columns, col_num)
    }
  }
  expr_pcs <- expr_pcs[, valid_columns]



  expr_pc_names <- colnames(expr_pcs)
  expr_pcs <- as.matrix(expr_pcs)
  num_expr_pcs <- dim(expr_pcs)[2]

    for (num_expr_pc in 1:num_expr_pcs) {
            expr_pc_vec <- expr_pcs[,num_expr_pc]
            expr_pc_name <- expr_pc_names[num_expr_pc]
            #print(paste0(num_pc, " - ", num_cov))
            lin_model <- lm(vec ~ expr_pc_vec)
            pve <- summary(lin_model)$adj.r.squared

            print(paste0(expr_pc_name, "\t", pve))
    }


}

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_covariate_pc_loading_correlation_heatmap <- function(expr_pcs, loadings, cluster_boolean) {
  loadings <- loadings[,2:(dim(loadings)[2])]
  expr_pcs <- expr_pcs[,2:(dim(expr_pcs)[2])]



  # Get stdev of each column of expr_pcs and remove columns with no variation
  valid_columns = c()
  for (col_num in 1:(dim(expr_pcs)[2])) {
    vary = var(expr_pcs[,col_num], na.rm=TRUE)
    if (vary != 0.0) {
      valid_columns <- c(valid_columns, col_num)
    }
  }
  expr_pcs <- expr_pcs[, valid_columns]



  expr_pc_names <- colnames(expr_pcs)

  loadings <- as.matrix(loadings)
  expr_pcs <- as.matrix(expr_pcs)



  # Initialize PVE heatmap
  factor_colnames <- paste0("PRS_PC", 1:(dim(loadings)[2]))
  factor_rownames <- expr_pc_names
  pve_map <- matrix(0, dim(expr_pcs)[2], dim(loadings)[2])

  colnames(pve_map) <- factor_colnames
  rownames(pve_map) <- factor_rownames


    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(loadings)[2]
    num_expr_pcs <- dim(expr_pcs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_expr_pc in 1:num_expr_pcs) {
            pc_vec <- loadings[,num_pc]
            expr_pc_vec <- expr_pcs[,num_expr_pc]
            #print(paste0(num_pc, " - ", num_cov))
            lin_model <- lm(pc_vec ~ expr_pc_vec)
            pve_map[num_expr_pc, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    


    if (cluster_boolean == TRUE) {
      ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

      melted_mat <- melt(pve_map)
      colnames(melted_mat) <- c("Expression_PC", "SURGE_PC","PVE")

      melted_mat$Expression_PC = factor(melted_mat$Expression_PC, levels=factor_rownames[ord])
      melted_mat$SURGE_PC = factor(melted_mat$SURGE_PC, levels=factor_colnames)
    } else {
      melted_mat <- melt(pve_map)
      colnames(melted_mat) <- c("Expression_PC", "SURGE_PC","PVE")

      melted_mat$Expression_PC = factor(melted_mat$Expression_PC, levels=factor_rownames)
      melted_mat$SURGE_PC = factor(melted_mat$SURGE_PC, levels=factor_colnames)      
    }

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Expression_PC, y=SURGE_PC)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + labs(y="",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)
}

scatter_plot_comparing_tcsc_h2_and_prs_weights <- function(df) {
    corry <- cor(df$prs_weights, df$tcsc_h2)
     plotter <- ggplot(df) + 
               geom_point(aes(x=prs_weights, y=tcsc_h2)) +
               gtex_v8_figure_theme() + 
               labs(x="PRS weights", y = "TCSC h2", title=paste0("spearman rho: ", corry)) + 
               theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) 

    return(plotter)
}


trait_name <- args[1]
bivariate_cafeh_output_dir <- args[2]
ukbb_prs_dir <- args[3]
analyzed_ukbb_prs_dir <- args[4]
output_stem <- args[5]



methods <- c("coloc", "causal_v1_coloc", "causal_v2_coloc", "mm_v1_coloc", "mm_v2_coloc")
for (method_iter in 1:length(methods)) {
method <- methods[method_iter]


p_active_thresholds <- c("0.1", "0.3", "0.5", "0.7", "0.9")

output_dir <- paste0(output_stem, trait_name, "_", method, "_")


####################################
# PRS weights at various thresholds
####################################
plot_arr <- list()
for (itera in 1:length(p_active_thresholds)) {
  p_active_threshold <- p_active_thresholds[itera]
  prs_weights_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_coloc_thresh_", p_active_threshold, "_", method, "_standardized_prs_weights_all_samples.txt")
  prs_weights <- read.table(prs_weights_file, header=TRUE, sep="\t")
  
  
  barplot <- make_prs_weight_bar_plot_with_standard_errors(prs_weights, p_active_threshold, "PRS weight")
  if (itera != length(p_active_thresholds)) {
      plot_arr[[itera]] <- barplot + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  } else {
    plot_arr[[itera]] <- barplot
  }
}
merged = plot_grid(plotlist=plot_arr, ncol=1, rel_heights=c(1, 1, 1,1,3.0))
output_file <- paste0(output_dir, "prs_weights.pdf")
ggsave(merged, file=output_file, width=7.2, height=8.0, units="in")

####################################
# Relative R-squared at various thresholds
####################################
plot_arr <- list()
for (itera in 1:length(p_active_thresholds)) {
  p_active_threshold <- p_active_thresholds[itera]
  r_squared_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_coloc_thresh_", p_active_threshold, "_", method, "_relative_r_squared.txt")
  r_squared_df <- read.table(r_squared_file, header=TRUE, sep="\t")
  
  
  barplot <- make_r_squared_bar_plot_with_standard_errors(r_squared_df, as.character(r_squared_df$prs_name))
  if (itera != length(p_active_thresholds)) {
      plot_arr[[itera]] <- barplot + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  } else {
    plot_arr[[itera]] <- barplot
  }
}
merged = plot_grid(plotlist=plot_arr, ncol=1, rel_heights=c(1, 1, 1,1,3.0))
output_file <- paste0(output_dir, "relative_r_squared_barplot.pdf")
ggsave(merged, file=output_file, width=7.2, height=8.0, units="in")



}






















if (FALSE) {



method = "adaptive_prior_coloc" 

output_dir <- paste0(output_dir, trait_name, "_", method, "_")
p_active_threshold <- "0.5"


####################################
# Adaptive prior coloc PRS weights at single threshold
####################################
prs_weights_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_coloc_thresh_", p_active_threshold, "_", method, "_standardized_prs_weights_all_samples.txt")
prs_weights <- read.table(prs_weights_file, header=TRUE, sep="\t")
prs_weights_barplot <- make_prs_weight_bar_plot_with_standard_errors(prs_weights, p_active_threshold, "PRS weight")

output_file <- paste0(output_dir, "prs_weights.pdf")
ggsave(prs_weights_barplot, file=output_file, width=7.2, height=4.0, units="in")

####################################
# Adaptive prior coloc PRS weights compared to TCSC (scatter plot)
####################################

tcsc_results_dir <- "/n/groups/price/tiffany/subpheno/AllGTExTissues_restore/1KG/STA_v8_320EUR_smalltoo/"
trait_tcsc_file <- paste0(tcsc_results_dir, "February_SingleTraitAnalysis_onejkweightperchunk_outliers_UKB_460K.", trait_name, "_2e+06_cislocus.txt.gz")
tcsc_df <- read.csv(trait_tcsc_file, sep="\t",header=T)

tissues <- as.character(prs_weights$prs_name)
tcsc_h2 <- c()
tcsc_h2_se <- c()

for (tissue_iter in 1:length(tissues)) {
  tissue_name = tissues[tissue_iter]
  tcsc_index = which(tcsc_df$Tissue == tissue_name)
  tcsc_h2 <- c(tcsc_h2, tcsc_df$Variance[tcsc_index])
  tcsc_h2_se <- c(tcsc_h2_se,tcsc_df$JK_SE[tcsc_index])
}

subset_tcsc_df <- data.frame(Tissue=tissues, tcsc_h2=tcsc_h2, tcsc_h2_se=tcsc_h2_se, prs_weights=prs_weights$weight)

scatter <- scatter_plot_comparing_tcsc_h2_and_prs_weights(subset_tcsc_df)
output_file <- paste0(output_dir, "prs_weights_tcsc_h2_scatter.pdf")
ggsave(scatter, file=output_file, width=7.2, height=4.0, units="in")

####################################
# Adaptive prior coloc PRS weights compared to TCSC (bar plot plot)
####################################
tcsc_barplot <- make_tcsc_h2_barplot_with_standard_errors(subset_tcsc_df)
merged = plot_grid(prs_weights_barplot + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()), tcsc_barplot, ncol=1, rel_heights=c(1,1.9))
output_file <- paste0(output_dir, "prs_weights_tcsc_h2_joint_barplot.pdf")

ggsave(merged, file=output_file, width=7.2, height=7.0, units="in")













}










































if (FALSE) {
###################################
# Num coloc components per tissue
###################################
num_components_per_tissue_file <- paste0(bivariate_cafeh_output_dir, method, "_results_", trait_name, "_0.5_num_prs_components.txt")
tissue_df = read.table(num_components_per_tissue_file, header=TRUE, sep="\t")
abs_barplot_5 <- make_num_coloc_component_bar_plot(tissue_df, paste0(trait_name , " 0.5"))
ratio_barplot_5 <- make_num_coloc_component_per_eqtl_component_bar_plot(tissue_df, paste0(trait_name , " 0.5"))

num_components_per_tissue_file <- paste0(bivariate_cafeh_output_dir, method, "_results_", trait_name, "_0.9_num_prs_components.txt")
tissue_df = read.table(num_components_per_tissue_file, header=TRUE, sep="\t")
abs_barplot_9 <- make_num_coloc_component_bar_plot(tissue_df, paste0(trait_name , " 0.9"))
ratio_barplot_9 <- make_num_coloc_component_per_eqtl_component_bar_plot(tissue_df, paste0(trait_name , " 0.9"))
joint <- plot_grid(abs_barplot_5, ratio_barplot_5, abs_barplot_9, ratio_barplot_9, ncol=2)
output_file <- paste0(output_dir, "joint_number_of_components.pdf")
ggsave(joint, file=output_file, width=7.2, height=8.0, units="in")


###################################
# Plot learned coloc priors
###################################
learned_coloc_prior_file <- paste0(bivariate_cafeh_output_dir, trait_name, "_learned_coloc_priors.txt")
coloc_prior_df <- read.table(learned_coloc_prior_file, header=FALSE, sep="\t", skip=1)
coloc_prior_barplot <- make_learned_coloc_prior_barplot(coloc_prior_df, trait_name)
output_file <- paste0(output_dir, "learned_coloc_prior_barplot.pdf")
ggsave(coloc_prior_barplot, file=output_file, width=7.2, height=4.0, units="in")


}

