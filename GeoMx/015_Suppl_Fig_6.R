# ## Script for Supplementary Figure 6 - DEGs within tumor heatmap
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(Biobase)
library(ComplexHeatmap)
library(circlize)

# Set the working directory
master.location <- setwd(master.location)

# Read dataset
Chord_tumor = readRDS('./analysis_files/5_Chord_tumor_clusters.RData')

# Load the data for the heatmap
load("./analysis_files/Suppl_Fig_6_Heatmap.RData")

#### Supplementary Figure 4 ####
# Plot the DEGs heatmap comparing tumor clusters
Heatmap(zscores[DEGs,],
        name = "Z-Score",
        border = TRUE,
        column_labels = pData(Chord_tumor)$ROI_ID,
        row_names_gp = gpar(fontsize = 6, fontface = "italic"),
        show_row_names = T, 
        show_column_names = F,
        row_split = names(DEGs),
        column_split = pData(Chord_tumor)$cluster,
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        column_title = "Segments",
        row_title = paste0("Top ", length(DEGs), " differentially expressed genes"),
        top_annotation = col_ann,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
