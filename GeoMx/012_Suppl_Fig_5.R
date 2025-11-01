## Script for Supplementary Figure 5 - Nuclei, Detected genes & T cells per cluster
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggpubr)

# Set working directory
master.location <- setwd(master.location)

# Load the data
load("./analysis_files/Suppl_Fig_5_boxplots.RData")

# Set the annotation colors
group_colors <- c("#1B9E77", "#D95F02")
names(group_colors) <- c("Inflamed", "Non-inflamed")

#### Supplementary Figure 5A ####
# Plot the number of captured nuclei per sample, for the stroma clusters
p1 <- ggplot(data = nuclei_stro, aes(x = groups, y = Nuclei, fill = groups))+
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 1500))+
  xlab("Stroma cluster")+
  ylab("Number of captured nuclei")+
  scale_x_discrete(limits = c("Non-inflamed", "Inflamed"))+
  scale_fill_manual(values = group_colors)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none")
p1

#### Supplementary Figure 5B ####
# Plot the number of detected genes per sample, for the stroma clusters
p2 <- ggplot(data = nuclei_stro, aes(x = groups, y = GenesDetected, fill = groups))+
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 7000))+
  xlab("Stroma cluster")+
  ylab("Number of detected genes")+
  scale_x_discrete(limits = c("Non-inflamed", "Inflamed"))+
  scale_fill_manual(values = group_colors)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none")
p2

#### Supplementary Figure 5C ####
# Plot the T cell count per sample, for the tumor clusters
p3 <- ggplot(data = IMC_counts, aes(x = groups_tumor, y = Tcells, fill = groups_tumor))+
  geom_boxplot()+
  xlab("Tumor cluster")+
  ylab("Mean T cell density / mm2")+
  scale_x_discrete(labels = c("Non-inflamed", "Inflamed"))+
  scale_fill_manual(values = group_colors) +
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
p3

#### Supplementary Figure 5D ####
# Plot the number of captured nuclei per sample, for the tumor clusters
p4 <- ggplot(data = nuclei_tumor, aes(x = groups, y = Nuclei, fill = groups))+
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 1500))+
  xlab("Tumor cluster")+
  ylab("Number of captured nuclei")+
  scale_x_discrete(limits = c("Non-inflamed", "Inflamed"))+
  scale_fill_manual(values = group_colors)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none")
p4

#### Supplementary Figure 5B ####
# Plot the number of detected genes per sample, for the tumor clusters
p5 <- ggplot(data = nuclei_tumor, aes(x = groups, y = GenesDetected, fill = groups))+
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 7000))+
  xlab("Tumor cluster")+
  ylab("Number of detected genes")+
  scale_x_discrete(limits = c("Non-inflamed", "Inflamed"))+
  scale_fill_manual(values = group_colors)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "none")
p5

#### Combine all plots ####
A <- ggarrange(p1, p4, p2, p5, nrow = 2, ncol = 2, align = "hv")
A
B <- ggarrange(A, 
               ggarrange(NULL, p3, NULL, ncol = 3, widths = c(1, 2, 1)),
               ncol = 1,
               heights = c(2, 1))
B
