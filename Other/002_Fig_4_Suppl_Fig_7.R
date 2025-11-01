## Script for Figure 4 and Supplementary Figure 7 - interferon-gamma treatment of chordoma cell lines
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Set working directory
master.location <- setwd(master.location)

# Read the data
# Normalized fold change
data <- read.delim("./input_files/IFNg_fold_changes.txt", header = TRUE)

# Ct values
Ct <- read.delim("./input_files/IFNg_housekeeping_genes.txt", header = TRUE)

# Set the colors for annotations
colors <- c("#E2D200", "#46ACC8", "#E58601", "#B40F20")

# Make factors for the visualization
data$Condition <- factor(data$Condition, levels = c("mock", "1dgIFN", "3dgIFN"))
data$Transcript <- factor(data$Transcript, levels = c("B2M", "TBXTall", "TBXT201/203", "TBXT202"))

Ct$Condition <- factor(Ct$Condition, levels = c("mock", "1dgIFN", "3dgIFN"))

# Average the Ct values of the genes of interest by using both housekeeping genes
data_avg <- data %>% group_by(Cell_line, Transcript, Condition) %>%
  summarise(avg_FC = mean(Fold_change), sd_FC = sd(Fold_change))

#### Figure 4 ####
# Plot the barplots for the averaged normalized fold change gene expression, for all conditions
ggplot(data = data_avg, aes(x = Condition, y = avg_FC, fill = Cell_line))+
  geom_errorbar(aes(ymin = avg_FC-sd_FC, ymax = avg_FC+sd_FC),
                position = position_dodge(0.9))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  ylab("Normalized fold change (2^-ddCt)")+
  xlab("IFNg treatment")+
  scale_x_discrete(labels = c("Mock", "1d", "3d"))+
  facet_wrap(~Transcript, nrow = 1)+
  scale_fill_manual(name = "Cell line", values = colors)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank())

#### Supplementary Figure 7 ####
# Plot the barplots for the housekeeping genes Ct values
ggplot(data = Ct, aes(x = Condition, y = Ct, fill = Cell_line))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  ylab("Average Ct value from triplicates")+
  xlab("IFNg treatment")+
  scale_x_discrete(labels = c("Mock", "1d", "3d"))+
  facet_wrap(~Control, nrow = 2)+
  scale_fill_manual(name = "Cell line", values = colors)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank())
