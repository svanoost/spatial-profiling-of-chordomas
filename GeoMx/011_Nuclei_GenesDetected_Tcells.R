## Script for extracting number of nuclei and detected genes per segment / and T cell counts from IMC data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(Biobase)

# Set working directory
master.location <- setwd(master.location)

# Load the data
Chord_tumor = readRDS("./analysis_files/5_Chord_tumor_clusters.RData")
Chord_stro = readRDS("./analysis_files/5_Chord_stroma_clusters.RData")

IMC_counts <- read.delim( "./input_files/IMC_counts.txt")

# Exclude L5436 - IMC data contains dedifferentiated areas. GeoMx data only conventional areas
IMC_counts <- IMC_counts[!IMC_counts$L_ID == "L5436",]
IMC_counts$Tcells <- as.numeric(IMC_counts$Tcells)
IMC_counts$groups_tumor <- factor(IMC_counts$groups_tumor, 
                                  levels = c("Non-inflamed", "Inflamed"))

#### Extract the phenotypic data - including nuclei and genes detected ####
# Stroma
nuclei_stro <- pData(Chord_stro)
nuclei_stro$cluster <- as.character(nuclei_stro$cluster)
nuclei_stro$groups <- "Inflamed"
nuclei_stro[nuclei_stro$cluster == "1", "groups"] <- "Non-inflamed"
# Tumor
nuclei_tumor <- pData(Chord_tumor)
nuclei_tumor$cluster <- as.character(nuclei_tumor$cluster)
nuclei_tumor$groups <- "Inflamed"
nuclei_tumor[nuclei_tumor$cluster == "1", "groups"] <- "Non-inflamed"

#### Statistics ####
t_tests <- data.frame(Segment = 0, Feature = 0, p_val = 0)
# Student's t-test to compare the number of genes detected per segment, per cluster
# Stroma
P <- t.test(nuclei_stro[nuclei_stro$groups == "Inflamed", "GenesDetected"], 
       nuclei_stro[nuclei_stro$groups == "Non-inflamed", "GenesDetected"])
# P = 7.288e-06
t_tests[1,] <- c("Stroma", "GenesDetected", P$p.value)

# Tumor
P <- t.test(nuclei_tumor[nuclei_tumor$groups == "Inflamed", "GenesDetected"], 
       nuclei_tumor[nuclei_tumor$groups == "Non-inflamed", "GenesDetected"])
# P = 0.02
t_tests[2,] <- c("Tumor", "GenesDetected", P$p.value)

# Student's t-test to compare the number of captured nuclei per segment, per cluster
# Stroma
P <- t.test(nuclei_stro[nuclei_stro$groups == "Inflamed", "Nuclei"], 
       nuclei_stro[nuclei_stro$groups == "Non-inflamed", "Nuclei"])
# P = 0.007
t_tests[3,] <- c("Stroma", "Nuclei", P$p.value)

# Tumor
P <- t.test(nuclei_tumor[nuclei_tumor$groups == "Inflamed", "Nuclei"], 
       nuclei_tumor[nuclei_tumor$groups == "Non-inflamed", "Nuclei"])
# P = 0.26
t_tests[4,] <- c("Tumor", "Nuclei", P$p.value)

# Student's t-test to compare the T cell counts per sample, per tumor cluster
P <- t.test(IMC_counts[IMC_counts$groups_tumor == "Inflamed", "Tcells"],
       IMC_counts[IMC_counts$groups_tumor == "Non-inflamed", "Tcells"])
# p = 0.016
t_tests[5,] <- c("Tumor", "Tcells", P$p.value)

#### Save the data ####
# t-test results
write.table(t_tests, file = "./output_files/Statistics_Nuclei_GenesDetected_Tcells.txt",
            row.names = FALSE, sep = "\t")

# Data for the boxplots - Supplementary Figure 5
save(nuclei_stro, nuclei_tumor, IMC_counts, file = "./analysis_files/Suppl_Fig_5_boxplots.RData")
