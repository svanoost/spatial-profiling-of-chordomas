# ## Script for Supplementary Figure 3 - Tumor vs Stroma DEG heatmap
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
Chord_dataset = readRDS('./analysis_files/4_Chord_dataset_NORM.RData')

# Load the data for the heatmap
load("./analysis_files/Suppl_Fig_3_Heatmap.RData")

#### Supplementary Figure 3 ####
# Plot the DEGs heatmap comparing tumor vs stroma
Heatmap(zscores[DEGs,],
        name = "Z-Score",
        border = T,
        column_labels = pData(Chord_dataset)$ROI_ID,
        row_names_gp = gpar(fontsize = 8, fontface = "italic"),
        show_row_names = T, 
        show_column_names = F,
        column_split = pData(Chord_dataset)$Segment,
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        column_title = "Segments",
        row_title = "Differentially expressed genes",
        top_annotation = col_ann,
        left_annotation = row_ann,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
