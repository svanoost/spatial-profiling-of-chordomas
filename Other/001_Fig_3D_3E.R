## Script for Figure 3D & 3E - PAS/dPAS stainings and Ki-67 IHC, comparing inflamed with non-inflamed samples
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# Set working directory
master.location <- setwd(master.location)

# Read the data
IHC <- read.delim("./input_files/IHC_PAS.txt", header = TRUE)

# Set the colors for annotations
cluster_colors <- c("#1B9E77", "#D95F02")
names(cluster_colors) <- c("Inflamed", "Non-inflamed")

PAS_colors <- c("#9932cc", "#ffaaff")
names(PAS_colors) <- c("High", "Low")

# Make a factor for the visualization
IHC$Tumor_factor <- factor(IHC$Tumor, levels = c("Non-inflamed", "Inflamed"))

# Calculate the percentage of PAS high and low samples per group
PAS <- IHC %>% group_by(Tumor_factor,PAS) %>%
  count() %>%
  group_by(Tumor_factor) %>%
  mutate(total = sum(n)) %>%
  mutate(perc = n/total*100)

#### Figure 3D ####
# Plot the stacked barplots for the PAS scores per group
A <- ggplot(data = PAS, aes(x = Tumor_factor, y = perc, fill = PAS))+
  geom_bar(stat = "identity", color = "black")+
  xlab("Tumor cluster")+
  ylab("Percentage (%)")+
  scale_fill_manual(name = "PAS", values = PAS_colors)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
A

#### Figure 3E ####
# Plot the boxplots for the Ki-67 positivity per 1000 tumor cells, per group
B <- ggplot(data = IHC, aes(x = Tumor_factor, y = KI67, fill = Tumor))+
  geom_boxplot()+
  xlab("Tumor cluster")+
  ylab("Ki-67 positivity (% of 1000 cells)")+
  scale_fill_manual(name = "Tumor cluster", values = cluster_colors)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
B

#### Combine the plots ####
fig3 <- ggarrange(A, B, nrow = 1,
                  widths = c(0.4, 0.5))
fig3

#### Statistics ####
# Student's t-test to compare the level of Ki-67 positivity per cluster
t.test(IHC[IHC$Tumor == "Inflamed", "KI67"],
       IHC[IHC$Tumor == "Non-inflamed", "KI67"])
# P = 0.048
