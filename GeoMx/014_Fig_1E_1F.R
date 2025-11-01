## Script for Figure 1E + 1F - Tumor clusters and DEGs
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Load the data
Chord_tumor = readRDS("./analysis_files/5_Chord_tumor_clusters.RData")
load("./analysis_files/Fig_1E_UMAPs.RData")
load("./analysis_files/Fig_1F_Volcanoplot.RData")

# Set the colors for annotation
indiv_colors = c("#A349A4", "#FFFF33", "#E7298A", "#091833", "#1B9E77", "#D95F02", "#7570B3",  "#66A61E", "#8DD3C7",
                 "#9F000F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
names(indiv_colors) = unique(pData(Chord_tumor)$Tags)

cluster_colors <- c("#1B9E77", "#D95F02")
names(cluster_colors) <- c("2", "1") # Inflamed and Non-inflamed

#### Figure 1E ####
# Plot the UMAPs
# Annotated for cluster ID
p1 <- ggplot(umap_df, aes(x = -UMAP1, y = UMAP2, color = snn)) +
  geom_point(size = 2) +
  scale_colour_manual(values = cluster_colors, 
                      labels = c("Non-inflamed", "Inflamed"),
                      name = "Cluster") +
  theme_bw()+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
p1
# Annotated for sample ID
p2 <- ggplot(umap_df, aes(x = -UMAP1, y = UMAP2, color = Individual)) +
  geom_point(size = 2) +
  scale_colour_manual(values = indiv_colors,
                      name = "Sample") +
  theme_bw()+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
p2

#### Figure 1F ####
# Plot the volcano plot comparing the tumor clusters
p3 <- ggplot(data = r_test, aes(x = -Estimate, y = -log10(`Pr(>|t|)`), color = result_category, 
                                label = Gene))+
  geom_vline(xintercept = c(-0.8, 0.8), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.005), linetype = "dashed")+
  geom_point()+
  geom_text_repel(data = subset(r_test, top_genes == "yes"),
                  point.padding = 0.15, color = "black",size=4,
                  min.segment.length = .1, box.padding = .2,
                  max.overlaps = 20, fontface = "italic")+
  xlab("Fold change")+
  ylab("Significance -log10(P)")+
  scale_color_manual(name = "Significance", values = c("#D95F02", "grey", "#1B9E77"), 
                     breaks = c("Non-inflamed", "NS or NE", "Inflamed"),
                     labels = c("Non-inflamed", "NS", "Inflamed"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")
p3

#### Combine the plots ####
fig1 <- ggarrange(ggarrange(p1, p2, nrow = 1),
                  p3, nrow = 2, heights = c(0.4, 0.6))
fig1
