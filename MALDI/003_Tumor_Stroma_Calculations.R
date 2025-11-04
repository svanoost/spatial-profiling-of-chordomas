## Script for processing the MALDI-MSI data into tumor and stroma regions
# R version 4.4.1

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(tidyr)
library(dplyr)
library(data.table)

# Set working directory
master.location <- setwd(master.location)

# Read the data
# Raw peak matrices
load("peakMat_merged.RData")

# UMAP and cluster metadata
UMAPimage <- fread("MALDI_UMAP_clustering_metadata.csv")

# load the identified lipid and metabolite names
peaks <- fread("./input_files/MALDI_peak_names.txt")$peaks

# Load the sample metadata with the GeoMx clusters
metadata <- fread("./input_files/MALDI_sample_metadata.txt")

#### Calculate tissue percentages ####
# Calculate the percentage of tumor and stroma spots for each sample
df <- UMAPimage[, c("Spot_ID", "L_ID", "tissue_type")] %>% 
  filter(tissue_type == "Tumor" | tissue_type == "Stroma") %>%
  group_by(L_ID) %>%
  count(tissue_type) %>%
  group_by(L_ID) %>%
  mutate(L_ID=L_ID, tissue_type=tissue_type, n=n, total_count=sum(n), perc=n/total_count*100)

# Pivot the percentages into separate columns based on tissue type
df_pivot <- pivot_wider(df[, c("L_ID", "tissue_type", "perc")], 
                        names_from = tissue_type, values_from = perc)
# Adjust the column names
colnames(df_pivot)[2:3] <- c("stroma_perc", "tumor_perc")

# Merge with the sample metadata
metadata <- left_join(metadata, df_pivot, by = "L_ID")

#### Calculate average TIC per tissue region per sample ####
# This is to compare tumor and stroma regions through differential analyses
# Tumor ####
# Add the tissue type to the peak matrices
tumor_calc <- left_join(peakMat, UMAPimage[,c("Spot_ID", "tissue_type")])

# Filter for tumor tissue
tumor_calc <- tumor_calc %>% filter(tissue_type == "Tumor")

# Filter for lipids and metabolites and group them per sample ID to calculate the average
summary_df <- tumor_calc %>% select(-Spot_ID, -x, -y, -xy) %>%
  group_by(L_ID) %>%
  summarise(across(starts_with("X"), mean, na.rm = TRUE))

# Stroma ####
# Add the tissue type to the peak matrices
stroma_calc <- left_join(peakMat, UMAPimage[,c("Spot_ID", "tissue_type")])

# Filter for stroma tissue
stroma_calc <- stroma_calc %>% filter(tissue_type == "Stroma")

# Filter for lipids and metabolites and group them per sample ID to calculate the average
summary_df_stroma <- stroma_calc %>% select(-Spot_ID, -x, -y, -xy) %>%
  group_by(L_ID) %>%
  summarise(across(starts_with("X"), mean, na.rm = TRUE))

# Set the column names of the averaged data sets according to the lipid and metabolite names
colnames(summary_df)[2:296] <- peaks
colnames(summary_df_stroma)[2:296] <- peaks

#### Save the data ####
# Tissue percentages -> integrate into sample metadata file
fwrite(df, file = "./input_files/MALDI_sample_metadata.txt", row.names = F, sep = "\t")

# Average TIC values for the tumor regions
fwrite(summary_df, file = "./input_files/MALDI_TumorIntensities.txt", row.names = F, sep = "\t")

# Average TIC values for the stroma regions
fwrite(summary_df_stroma, file = "./input_files/MALDI_StromaIntensities.txt", , row.names = F, sep = "\t")
