# ## Script for Figure S6A + S6B - DEPs within stroma and tumor compartments between inflammatory clusters
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ComplexHeatmap)
library(circlize)

# Set the working directory
master.location <- setwd(master.location)

# Read the GeoMx data per compartment
Chord_stro = readRDS('./analysis_files/5_Chord_stroma_clusters.RData')
Chord_tum = readRDS('./analysis_files/5_Chord_tumor_clusters.RData')

# Load the data for the heatmaps
load("./analysis_files/Fig_S6A_Heatmap.RData")
load("./analysis_files/Fig_S6B_Heatmap.RData")

#### Figure S6A ####
# Generate the DEPs heatmap comparing the stroma clusters
h1 <- Heatmap(zscores_ssGSEA[DEPs,],
        name = "Z-Score",
        border = T,
        row_names_gp = gpar(fontsize = 8),
        show_row_names = T, 
        show_column_names = F,
        column_split = pData(Chord_stro)$cluster,
        row_split = names(DEPs),
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        column_title = "Segments",
        row_title = paste0("Top 10% differentially enriched pathways"),
        top_annotation = col_ann,
        heatmap_legend_param = list(
          direction = "horizontal"
        ),
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

# Plot the heatmap
draw(h1, annotation_legend_side = "bottom", heatmap_legend_side = "bottom", 
     merge_legend = TRUE)

#### Figure S6B ####
# Generate the DEPs heatmap comparing the tumor clusters
h2 <- Heatmap(zscores_ssGSEA_tumor[DEPs_tumor,],
              name = "Z-Score",
              border = T,
              row_names_gp = gpar(fontsize = 8),
              show_row_names = T, 
              show_column_names = F,
              column_split = pData(Chord_tum)$cluster,
              row_split = names(DEPs_tumor),
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_title = "Segments",
              row_title = paste0("Top 10% differentially enriched pathways"),
              top_annotation = col_ann_tumor,
              heatmap_legend_param = list(
                direction = "horizontal"
              ),
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

# Plot the heatmap
draw(h2, annotation_legend_side = "bottom", heatmap_legend_side = "bottom", 
     merge_legend = TRUE)
