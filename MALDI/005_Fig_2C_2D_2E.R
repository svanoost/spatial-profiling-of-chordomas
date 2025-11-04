## Script for processing the MALDI-MSI data into tumor and stroma regions
# R version 4.4.1

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Read the data
# Tumor vs Stroma regions
load("./analysis_files/MALDI_Tumor_vs_Stroma_regions.RData")

# Inflamed vs Non-inflamed regions
load("./analysis_files/MALDI_Inflamed_vs_NonInflamed_regions.RData")

#### Figure 2C ####
# Plot the volcano plot comparing tumor with stroma regions
C <- ggplot(data = DF, aes(x = LogFC, y = -log10(Padj), color = diffexpressed, 
                           label = delabel))+
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_point()+
  geom_text_repel(data = subset(DF[!grepl("X", DF$delabel),], diffexpressed != "Not significant"),
                  point.padding = 0.15, color = "black",size = 3,
                  min.segment.length = .1, box.padding = .2,
                  max.overlaps = 20)+
  xlab("Log2 fold change")+ ylab("Significance -log10(P)")+
  scale_color_manual(name = "Significance", values = c('#663300', "grey", '#f17a7b'), 
                     breaks = c("Downregulated", "Not significant", "Upregulated"),
                     labels = c("Stroma", "NS", "Tumor"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")
C

#### Figure 2D ####
# Plot the volcano plot comparing inflamed with non-inflamed stroma regions
D <- ggplot(data = Stroma_df, aes(x = StromaLogFC, y = -log10(StromaPval), color = diffexpressed, 
                                  label = delabel))+
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_point()+
  geom_text_repel(data = Stroma_df[!grepl("X", Stroma_df$delabel) & Stroma_df$StromaPval < 0.05,],
                  point.padding = 0.15, color = "black",size = 3,
                  min.segment.length = .1, box.padding = .2,
                  max.overlaps = 20)+
  xlab("Log2 fold change")+ ylab("Significance -log10(P)")+
  scale_color_manual(name = "Significance", values = c('#1B9E77', "grey", '#D95F02'), 
                     breaks = c("Upregulated", "Not Significant", "Downregulated"),
                     labels = c("Inflamed", "NS", "Non-Inflamed"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")
D

#### Figure 2E ####
# Plot the volcano plot comparing inflamed with non-inflamed tumor regions
E <- ggplot(data = Tumor_df, aes(x = TumorLogFC, y = -log10(TumorPval), color = diffexpressed, 
                                 label = delabel))+
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_point()+
  geom_text_repel(data = Tumor_df[!grepl("X", Tumor_df$delabel) & Tumor_df$TumorPval < 0.05,],
                  point.padding = 0.15, color = "black",size = 3,
                  min.segment.length = .1, box.padding = .2,
                  max.overlaps = 20)+
  xlab("Log2 fold change")+ylab("Significance -log10(P)")+
  scale_color_manual(name = "Significance", values = c('#1B9E77', "grey", '#D95F02'), 
                     breaks = c("Upregulated", "Not Significant", "Downregulated"),
                     labels = c("Inflamed", "NS", "Non-Inflamed"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")
E

#### Combine the plots ####
# First combine the inflamed vs non-inflamed plots
D_E <- ggarrange(D, E, nrow = 1)

# Then add the Tumor vs Stroma plot
combined <- ggarrange(ggarrange(C, NULL, nrow = 1), 
                      D_E,
                  ncol = 1, heights = c(0.6, 0.4))
combined
