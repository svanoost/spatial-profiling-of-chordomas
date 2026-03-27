## Script for Figure 2A and 2B - UMAP and clusters of the MALDI-MSI data and their spatial localization
# R version 4.4.1

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(data.table)
library(dplyr)
library(cowplot)

# Set working directory
master.location <- setwd(master.location)

# Read the data
UMAP <- fread("./input_files/MALDI_UMAP_clustering_metadata.txt")
UMAP$clusters <- as.character(UMAP$clusters)

# Set the colors for annotation
# Samples from the same patient are attributed the same color in UMAP_colors[["Patient"]]
UMAP_colors <- list(clusters = c("0" = "#663300", "1" = "#005300", "2" = "#009FFF", "3" = "#FF00B6", "4" = "#00FF00",
                                "5" = "#0000FF", "6" = "#FFD300", "7" = "#FF0000", "8" = "#000033",
                                "9" = "#783FC1", "10" = "#00FFBE"),
                   Patient = c("3034" = "#BEBADA", "3356" = "#FFFF33", "3382" = "#A349A4", "3519" = "#1B9E77",
                               "3770" = "#091833", "4716" = "#FDB462", "5327" = "#D95F02", "5343" = "#7570B3",
                               "5436" = "#E7298A", "6152" = "#8DD3C7", "6325" = "#80B1D3", "6371" = "#9F000F",
                               "6473" = "#E6AB02", "6737" = "#FB8072", "6951" = "#B3DE69", "985" = "#66A61E",
                               "6448" = "#7570B3", "2038" = "#66A61E", "4336" = "#66A61E", "6986" = "#FCCDE5"))

# Create a new column to extract the sample ID (which is actual L_ID: e.g. L3034 etc.)
UMAP$Patient <- UMAP$L_ID

# Convert to data.table if it's not already
setDT(UMAP)
# Split the L_ID column based on "_" and keep the 3rd element
UMAP[, Patient := tstrsplit(L_ID, "_", keep = 3)]
# Check that indeed the sample IDs match
unique(UMAP$Patient)

# Create a factor for the visualization of the legend
UMAP$clusters <- factor(UMAP$clusters, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

#### Figure 2A ####
# Plot the UMAP and MALDI-MSI clusters
# Use rasterise to decrease the size of the image. Image was saved with pdf()
#pdf(file = "UMAP_clusters_raster_w50_h50.pdf, width = 50, height = 50)
ggplot(data = UMAP, aes(x = UMAP1, y = UMAP2, color = clusters))+
  rasterise(geom_point(size = 0.2))+
  scale_color_manual(values = UMAP_colors$clusters)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none")
#dev.off()

# Due to the small point size of the data points in the UMAP, the legend was plotted separately
plot <- ggplot(data = UMAP, aes(x = UMAP1, y = UMAP2, color = clusters))+
  rasterise(geom_point(size = 0.2))+
  scale_color_manual(values = UMAP_colors$clusters, name = "Clusters")+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank())+
  guides(
  color = guide_legend(
    override.aes = list(size = 2),
    ncol = 2
  ))

# Extract the legend and plot it as a ggplot
legend_umap <- get_legend(plot)
as_ggplot(legend_umap)

#### Figure 2B ####
# Plot the spatial localization of the clusters - to match with the corresponding H&E images (from the same slides)
# Samples L4336, L6448, and L6371 were selected for visualization
to_vis <- c("6371", "6448", "4336")
# Use rasterise to decrease the size of the image. Image was saved with pdf()
#pdf(file = "Spatial_localization_of_clusters_examples.pdf", width = 5, height = 15)
ggplot(data = UMAP[UMAP$Patient %in% to_vis,], aes(x = x, y = -y, fill = clusters))+
  rasterize(geom_tile(), dpi = 300)+
  scale_fill_manual(values = UMAP_colors$clusters)+
  facet_wrap(~Patient, ncol = 1)+
  theme_void()+
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
#dev.off()
