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
               geom_point(aes(x=num_eqtl_components, y=num_cafeh_components)) +
               gtex_v8_figure_theme() + 
               labs(x="Number of eQTL CAFEH components", y = "Number of colocalizing CAFEH components") + 
               theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
    return(plotter)
}

make_tissue_sample_size_num_cafeh_component_scatter <- function(tissue_df) {
    plotter <- ggplot(tissue_df) + 
               geom_point(aes(x=sample_size, y=num_cafeh_components)) +
               gtex_v8_figure_theme() + 
               labs(x="Tissue sample size", y = "Number of CAFEH components") + 
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

make_prs_weight_bar_plot_with_standard_errors <- function(prs_weights, ordered_studies) {
  prs_weights$prs_name = factor(prs_weights$prs_name, levels=ordered_studies)
    p <- ggplot(data=prs_weights, aes(x=prs_name, y=weight)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(legend.position="top") +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    labs(y="PRS NNLS weight", x="", fill="") +
    geom_errorbar(aes(ymin=weight-(weight_standard_error), ymax=weight+(weight_standard_error)), position = position_dodge(), width = .75, size=.2)

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
  tissue_df$fraction = tissue_df$num_cafeh_components/tissue_df$num_eqtl_components
  print(head(tissue_df))
  p <- ggplot(data=tissue_df, aes(x=tissue, y=fraction)) +
    geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
    gtex_v8_figure_theme() +
    theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    labs(y="Fraction of eQTL components that coloc", x="", fill="")
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




trait_name <- args[1]
bivariate_cafeh_output_dir <- args[2]
ukbb_prs_dir <- args[3]
analyzed_ukbb_prs_dir <- args[4]
output_dir <- args[5]
model_version <- args[6]

output_dir <- paste0(output_dir, trait_name, "_", model_version, "_")

num_components_per_tissue_file <- paste0(bivariate_cafeh_output_dir, "cafeh_results_", trait_name, "_all_tissues_0.5_num_prs_components.txt")

# Load in PRS scores in each tissue
tissue_specific_prs_scores_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_standardized_residual_prs_scores.txt")
tissue_specific_prs_scores <- read.table(tissue_specific_prs_scores_file, header=TRUE, sep="\t")
scores_mat = as.matrix(tissue_specific_prs_scores[,2:(dim(tissue_specific_prs_scores)[2])])

# Load in MOR posteriors
#mor_posterior_file <- paste0(input_dir, trait_name, "_mixture_of_regressions_posteriors.txt")
#mor_posterior <- read.table(mor_posterior_file, header=TRUE, sep="\t")


# Load in learned weights for each tissue prs
prs_weights_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_prs_weights_all_samples.txt")
prs_weights <- read.table(prs_weights_file, header=TRUE, sep="\t")


# Compute global prs
global_prs = (scores_mat %*% prs_weights$weight)[,1]

# Loaded in R-squared results for each tissue
relative_r_squared_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_relative_r_squared.txt")
relative_r_squared <- read.table(relative_r_squared_file, header=TRUE, sep="\t")

tissue_relative_r_squared <- relative_r_squared[1:23,]


# Loaded in R-squared results for each tissue
relative_pca_r_squared_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_relative_r_squared_with_prs_pcs.txt")
relative_pca_r_squared <- read.table(relative_pca_r_squared_file, header=TRUE, sep="\t")


# Correlation of PRS weights between studies
corr_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_prs_correlation_matrix.txt")
corr_mat <- read.table(corr_file, header=TRUE, sep="\t")
corr_mat <- corr_mat[,2:(dim(corr_mat)[2])]

# Pca loadings
prs_pca_loadings_file = paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_prs_pca_loadings.txt")
prs_pca_loadings = read.table(prs_pca_loadings_file, header=TRUE, sep="\t")

# Pca components
prs_pca_components_file = paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_prs_pca_principal_components.txt")
prs_pca_pcs = read.table(prs_pca_components_file, header=TRUE, sep="\t")

# PCA VE 
prs_pca_ve_file = paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version,"_prs_pca_variance_explained.txt")
prs_pca_ve = read.table(prs_pca_ve_file, header=TRUE, sep="\t")

# Get number of components per tissue
tissue_df = read.table(num_components_per_tissue_file, header=TRUE, sep="\t")
print(head(tissue_df))
sample_size_order = order(-tissue_df$sample_size)
num_comp_order = order(-tissue_df$num_cafeh_components)
num_eqtl_order = order(-tissue_df$num_eqtl_components)
tissues_ordered_alphabetically = as.character(tissue_df$tissue)
tissues_ordered_by_sample_size = as.character(tissue_df$tissue[sample_size_order])
tissues_ordered_by_number_cafeh_components = as.character(tissue_df$tissue[num_comp_order])
tissues_ordered_by_number_eqtl_components = as.character(tissue_df$tissue[num_eqtl_order])


# Prs weight/num eqtl
prs_weights_per_eqtl_comp = data.frame(prs_weights)
prs_weights_per_eqtl_comp$weight = prs_weights_per_eqtl_comp$weight/tissue_df$num_eqtl_components
prs_weights_per_eqtl_comp$weight_standard_error = prs_weights_per_eqtl_comp$weight_standard_error/tissue_df$num_eqtl_components

# Relative R-squared/num eqtl
tissue_relative_r_squared_per_eqtl_comp = data.frame(tissue_relative_r_squared)
tissue_relative_r_squared_per_eqtl_comp$relative_r_squared = tissue_relative_r_squared_per_eqtl_comp$relative_r_squared/tissue_df$num_eqtl_components
tissue_relative_r_squared_per_eqtl_comp$relative_r_squared_standard_error = tissue_relative_r_squared_per_eqtl_comp$relative_r_squared_standard_error/tissue_df$num_eqtl_components

# Blood-related covariates
blood_cov_file = paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_blood_covariates.txt")
blood_cov = read.table(blood_cov_file, header=TRUE, sep="\t")

# Trait related covariates
trait_cov_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_trait_covariates.txt")
trait_cov = read.table(trait_cov_file, header=TRUE, sep="\t")

# Technical covariates
technical_cov_file <- paste0(analyzed_ukbb_prs_dir, trait_name, "_", model_version, "_technical_covariates.txt")
technical_cov = read.table(technical_cov_file, header=TRUE, sep="\t")



#temp_correlation(mor_posterior$mor_component_1, blood_cov)

###################
# Scatter plot correlating gtex tissue sample size with number of cafeh components identified
###################
output_file <- paste0(output_dir, "tissue_sample_size_vs_num_cafeh_components_scatter.pdf")
scatter <- make_tissue_sample_size_num_cafeh_component_scatter(tissue_df)
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")


###################
# Scatter plot correlating num cafeh eqtl components with number of cafeh components identified
###################
output_file <- paste0(output_dir, "cafeh_eqtl_components_vs_cafeh_coloc_components_scatter.pdf")
scatter <- make_cafeh_eqtl_component_vs_cafeh_coloc_component_scatter(tissue_df)
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")



###################
# Bar plot showing fraction of eqtl components that colocalize
###################
output_file <- paste0(output_dir, "fraction_of_eqtl_components_that_colocalize.pdf")
bar <- make_fraction_of_eqtl_components_that_colocalize(tissue_df)
ggsave(bar, file=output_file, width=7.2, height=6.0, units="in")


###################
# Heatmap showing correlation of PRS scores between tissues
###################
heatmap <- make_heatmap(corr_mat)
ggsave(heatmap, file=paste0(output_dir, "tissue_correlation_heatmap.pdf"), width=7.2, height=6.0, units="in")

###################
# Histogram showing correlation of PRS scores between pairs of tissues
###################
histo <- make_pairwise_prs_correlation_histogram(corr_mat)
ggsave(histo, file=paste0(output_dir, "tissue_pairwise_correlation_histogram.pdf"), width=7.2, height=4.0, units="in")


###################
# Bar plot showing r-squared of each of prs models
###################
output_file <- paste0(output_dir, "relative_r_squared_bar_plot_ordered_by_num_eqtl_components.pdf")
ordered_studies <- c("joint_prs", tissues_ordered_by_number_eqtl_components)
barplot <- make_r_squared_bar_plot_with_standard_errors(relative_r_squared, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "relative_r_squared_bar_plot_ordered_by_number_of_cafeh_components.pdf")
ordered_studies <- c("joint_prs", tissues_ordered_by_number_cafeh_components)
barplot <- make_r_squared_bar_plot_with_standard_errors(relative_r_squared, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "relative_r_squared_bar_plot_ordered_alphabetically.pdf")
ordered_studies <- c("joint_prs", tissues_ordered_alphabetically)
barplot <- make_r_squared_bar_plot_with_standard_errors(relative_r_squared, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "relative_r_squared_bar_plot_ordered_alphabetically_tissue_only.pdf")
ordered_studies <- tissues_ordered_alphabetically
barplot <- make_r_squared_bar_plot_with_standard_errors(relative_r_squared, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")



output_file <- paste0(output_dir, "relative_r_squared_bar_plot_ordered_by_num_eqtl_components_tissue_only.pdf")
ordered_studies <- tissues_ordered_by_number_eqtl_components
barplot <- make_r_squared_bar_plot_with_standard_errors(tissue_relative_r_squared, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "relative_r_squared_per_eqtl_component_bar_plot_ordered_by_num_eqtl_components_tissue_only.pdf")
ordered_studies <- tissues_ordered_by_number_eqtl_components
barplot <- make_r_squared_bar_plot_with_standard_errors(tissue_relative_r_squared_per_eqtl_comp, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "relative_r_squared_per_eqtl_component_bar_plot_ordered_alphabetically_tissue_only.pdf")
ordered_studies <- tissues_ordered_alphabetically
barplot <- make_r_squared_bar_plot_with_standard_errors(tissue_relative_r_squared_per_eqtl_comp, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")


###################
# Bar plot showing r-squared of each of PCA-prs models
###################
output_file <- paste0(output_dir, "relative_r_squared_pca_bar_plot.pdf")
barplot <- make_r_squared_bar_plot_with_standard_errors(relative_pca_r_squared, as.character(relative_pca_r_squared$prs_name))
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")



###################
# Bar plot showing prs weight bar plot
###################
output_file <- paste0(output_dir, "prs_weight_bar_plot_ordered_by_number_of_eqtl_components.pdf")
ordered_studies <- c(tissues_ordered_by_number_eqtl_components)
barplot <- make_prs_weight_bar_plot_with_standard_errors(prs_weights, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_weight_bar_plot_ordered_by_number_of_cafeh_components.pdf")
ordered_studies <- c(tissues_ordered_by_number_cafeh_components)
barplot <- make_prs_weight_bar_plot_with_standard_errors(prs_weights, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_weight_bar_plot_ordered_alphabetically.pdf")
ordered_studies <- c(tissues_ordered_alphabetically)
barplot <- make_prs_weight_bar_plot_with_standard_errors(prs_weights, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")

###################
# Bar plot showing prs weight per eqtl comp bar plot
###################
output_file <- paste0(output_dir, "prs_weight_per_eqtl_component_bar_plot_ordered_by_number_of_eqtl_components.pdf")
ordered_studies <- c(tissues_ordered_by_number_eqtl_components)
barplot <- make_prs_weight_per_eqtl_component_bar_plot_with_standard_errors(prs_weights_per_eqtl_comp, ordered_studies)
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")



if (FALSE) {
tcsc_h2 <- data.frame(prs_weights)
print(tcsc_h2)
tcsc_h2$weight <- c(.019, -.004, .044, -.006, .0010, .019, .0156, .002, -.006, .034, -.008, -.005)
tcsc_h2$weight_standard_error <- c(.008, .006, .012, .011, .010, .019, .008, .016, .01, .02, .008,.011)
output_file <- paste0(output_dir, "tcsc_h2_bar_plot_ordered_by_number_of_cafeh_components.pdf")
ordered_studies <- c(tissues_ordered_by_number_cafeh_components)
barplot <- make_prs_weight_bar_plot_with_standard_errors(tcsc_h2, ordered_studies) + labs(y="TCSC h2", x="", fill="") 
ggsave(barplot, file=output_file, width=7.2, height=6.0, units="in")
}

###################
# PRS PCA principal components heatmap
###################
output_file <- paste0(output_dir, "prs_pca_principal_component_heatmap.pdf")
heatmap <- make_prs_pca_principal_component_heatmap(prs_pca_pcs, ordered_studies)
ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")

###################
# Scatter plot showing correlation of PC1 and number of cafeh components
###################
output_file <- paste0(output_dir, "prs_pc1_vs_tissue_cafeh_components_scatter.pdf")
scatter <- prs_principal_component_number_of_components_scatter(prs_pca_pcs$prs_pc1, prs_pca_pcs$sample_name, tissue_df, 1)
ggsave(scatter, file=output_file, width=7.2, height=4.5, units="in")


#######################################
# PVE plot showing fraction of variance explained through each factor
#######################################
output_file <- paste0(output_dir, "PRS_pca_variance_explained_line_plot.pdf")
pve_plot <- make_pc_variance_explained_line_plot(prs_pca_ve)
ggsave(pve_plot, file=output_file, width=7.2, height=3.5, units="in")


###################
# Correlation heatmap betweeen PRS PCs and covariates
###################
# Blood covariates
output_file <- paste0(output_dir, "prs_pc_blood_covariate_correlation_heatmap.pdf")
#heatmap <- make_covariate_pc_loading_correlation_heatmap(blood_cov, prs_pca_loadings, TRUE)
#ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")

# Trait covariates
output_file <- paste0(output_dir, "prs_pc_trait_covariate_correlation_heatmap.pdf")
#heatmap <- make_covariate_pc_loading_correlation_heatmap(trait_cov, prs_pca_loadings, TRUE)
#ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")

# Technical covariates
output_file <- paste0(output_dir, "prs_pc_technical_covariate_correlation_heatmap.pdf")
#heatmap <- make_covariate_pc_loading_correlation_heatmap(technical_cov, prs_pca_loadings, FALSE)
#ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")


output_file <- paste0(output_dir, "whole_blood_prs_vs_brain_cortex_prs_colored_by_global_prs.pdf")
scatter <- make_prs_pc_scatter_colored_by_covariate(tissue_specific_prs_scores$Whole_Blood, tissue_specific_prs_scores$Brain_Cortex, global_prs, "Whole Blood PRS", "Brain Cortex PRS", "Global PRS")
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

#output_file <- paste0(output_dir, "whole_blood_prs_vs_brain_cortex_prs_colored_by_mor_posterior.pdf")
#scatter <- make_prs_pc_scatter_colored_by_covariate(tissue_specific_prs_scores$Whole_Blood, tissue_specific_prs_scores$Brain_Cortex, mor_posterior$mor_component_1, "Whole Blood PRS", "Brain Cortex PRS", "MOR posterior 1")
#ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

#output_file <- paste0(output_dir, "whole_blood_prs_vs_adipose_subcutaneous_prs_colored_by_mor_posterior.pdf")
#scatter <- make_prs_pc_scatter_colored_by_covariate(tissue_specific_prs_scores$Whole_Blood, tissue_specific_prs_scores$Adipose_Subcutaneous, mor_posterior$mor_component_1, "Whole Blood PRS", "Adipose Subcutaneous PRS", "MOR posterior 1")
#ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "whole_blood_prs_vs_adipose_subcutaneous_prs_density_plot.pdf")
scatter <- make_scatter_density_plot(tissue_specific_prs_scores$Whole_Blood, tissue_specific_prs_scores$Adipose_Subcutaneous, "Whole Blood PRS", "Adipose_Subcutaneous PRS")
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "whole_blood_prs_vs_brain_cortex_prs_density_plot.pdf")
scatter <- make_scatter_density_plot(tissue_specific_prs_scores$Whole_Blood, tissue_specific_prs_scores$Brain_Cortex, "Whole Blood PRS", "Brain Cortex PRS")
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")













# Make PRS PC scatter colored by covariate
output_file <- paste0(output_dir, "prs_pc1_vs_prs_pc2_scatter_colored_by_blood_white_count.pdf")
#scatter <- make_prs_pc_scatter_colored_by_covariate(prs_pcs$prs_pc1, prs_pcs$prs_pc2, cov$blood_WHITE_COUNT, "PRS_PC1", "PRS_PC2", "blood_white_count")
#ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_pc2_vs_prs_pc3_scatter_colored_by_blood_white_count.pdf")
#scatter <- make_prs_pc_scatter_colored_by_covariate(prs_pcs$prs_pc2, prs_pcs$prs_pc3, cov$blood_WHITE_COUNT, "PRS_PC2", "PRS_PC3", "blood_white_count")
#ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_pc1_vs_blood_white_count_scatter_colored_by_prs_pc2.pdf")
#scatter <- make_prs_pc_scatter_colored_by_covariate(prs_pcs$prs_pc1, cov$blood_WHITE_COUNT, prs_pcs$prs_pc2, "PRS_PC1", "blood_white_count", "PRS_PC2")
#ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_pc2_vs_blood_NEUTROPHIL_PCT.pdf")
#scatter <- make_scatter_density_plot(prs_pcs$prs_pc2, cov$blood_NEUTROPHIL_PCT, "PRS_PC2", "blood_neutrophil_pct")
#ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(output_dir, "prs_pc1_vs_blood_WHITE_COUNT.pdf")
#scatter <- make_scatter_density_plot(prs_pcs$prs_pc1, cov$blood_WHITE_COUNT, "PRS_PC1", "blood_white_count")
#ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")


#heatmap <- make_heatmap(corr_mat)
#ggsave(heatmap, file=paste0(output_dir, "headmap.pdf"))

