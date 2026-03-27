## Script for Figure 1A + 1B - Tumor vs Stroma DEGs and DEPs
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Set working directory
master.location <- setwd(master.location)

# Load the data
load("./analysis_files/Fig_1A_Volcanoplot.RData")
load("./analysis_files/Fig_1B_barplot.RData")

#### Figure 1A ####
# Plot the volcano plot - DEGs Tumor vs Stroma
p1 <- ggplot(data = r_test, aes(x = -Estimate, y = -log10(`Pr(>|t|)`), color = result_category, 
                                label = Gene))+
  geom_vline(xintercept = c(-0.7, 0.7), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.0001), linetype = "dashed")+
  geom_point(size = 2)+
  geom_text_repel(data = subset(r_test, top_genes == "yes"),
                  point.padding = 0.15, color = "black",size=3,
                  min.segment.length = .1, box.padding = .2,
                  max.overlaps = 20, fontface = "italic")+
  xlab("Fold change")+ylab("Significance -log10(P)")+
  scale_color_manual(name = "Significance", values = c('#663300', "grey", '#f17a7b'), 
                     breaks = c("Stroma", "NS or NE", "Tumor"),
                     labels = c("Stroma", "NS", "Tumor"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")
p1

#### Figure 1B ####
# Plot the bar plot - DEPs Tumor vs Stroma
p2 <- ggplot(data = top_paths, aes(x = -Estimate, y = pathway_factor, fill = result_category))+
  geom_bar(stat = "identity")+
  xlab("NES")+
  ylab("Top 20 enriched pathways")+
  scale_fill_manual(name = "Segment", breaks = c("Stroma", "Tumor"),
                    values = c("#663300", "#f17a7b"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none")
p2

#### Combine both plots ####
fig1 <- ggarrange(p1, p2, nrow = 1, widths = c(0.4, 0.6))
fig1
