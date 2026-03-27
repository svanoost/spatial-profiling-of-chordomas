## Script for the differential analyses with the MALDI-MSI data
# R version 4.4.1

#### Set up environment ####
rm(list = ls())

# Set working directory
master.location <- setwd(master.location)

# Read the data
# Average Tumor TIC values
Tumor <- read.delim("./input_files/MALDI_TumorIntensities.txt", header = T)[, 2:296]

# Average Stroma TIC values
Stroma <- read.delim("./input_files/MALDI_StromaIntensities.txt", header = T)[, 2:296]

# Sample metadata including the GeoMx clusters and the tissue percentages
metadata <- read.delim("./input_files/MALDI_sample_metadata.txt", header = T)

# Set the row names of the data sets according to the sample IDs
row.names(Tumor) <- metadata$L_ID
row.names(Stroma) <- metadata$L_ID

# Exclude L6473, since this sample was not included in the GeoMx analysis
metadata <- metadata[metadata$L_ID != "peakmat_ITO1_6473_T1",]
Tumor <- Tumor[row.names(Tumor) != "peakmat_ITO1_6473_T1",]
Stroma <- Stroma[row.names(Stroma) != "peakmat_ITO1_6473_T1",]

#### Differential analyses - Inflamed vs Non-inflamed ####
# Filter the GeoMx samples
# Tumor regions
TumorAverage <- Tumor[metadata$Stroma %in% c("Inflamed", "Non-inflamed"),]
# 6 Inflamed vs 7 Non-inflamed

# Stroma regions
StromaAverage <- Stroma[metadata$Stroma %in% c("Inflamed", "Non-inflamed"),]
# 6 Inflamed vs 7 Non-inflamed
# But, some samples have very little stroma, therefore, only take along samples > 2% of stroma in the sample
Inflam <- metadata[metadata$Stroma %in% c("Inflamed", "Non-inflamed"),]
Inflam_StromaFilter <- Inflam[Inflam$stroma_perc > 2,]
StromaAverage <- StromaAverage[Inflam$stroma_perc > 2,]
# 6 Inflamed vs 6 Non-inflamed

# Initialize for the loop
TumorPval <- NA
StromaPval <- NA
TumorFC <- NA
StromaFC <- NA

# Perform the t-tests and calculate P values and the fold changes (FC) for each lipid and metabolite
for(i in 1:295){
  # First, t-tests and P values for
  # Tumor regions
  TumorPval[i] <- t.test(TumorAverage[Inflam$Stroma == "Inflamed", i],
                       TumorAverage[Inflam$Stroma == "Non-inflamed", i])$p.val
  # Stroma regions
  StromaPval[i] <- t.test(StromaAverage[Inflam_StromaFilter$Stroma == "Inflamed", i],
                        StromaAverage[Inflam_StromaFilter$Stroma == "Non-inflamed", i])$p.val
  # Second, FCs for
  # Tumor regions
  TumorFC[i] <- mean(TumorAverage[Inflam$Stroma == "Inflamed", i]) /
    mean(TumorAverage[Inflam$Stroma == "Non-inflamed", i])
  # Stroma regions
  StromaFC[i] <- mean(StromaAverage[Inflam_StromaFilter$Stroma == "Inflamed", i]) /
    mean(StromaAverage[Inflam_StromaFilter$Stroma == "Non-inflamed", i])
}

# Calculate the False Discovery Rate (FDR) using the Benjamini-Hochberg method
TumorPadj <- p.adjust(TumorPval, method = "BH")
StromaPadj <- p.adjust(StromaPval, method = "BH")

# Create the data frames for visualization of the results
# Tumor regions ####
Tumor_df <- data.frame(TumorFC, TumorPval, TumorPadj)

# Make a column for differential abudance
Tumor_df$diffexpressed <- "Not significant"

# Log2 transform the FC
TumorLogFC <- log(TumorFC, base = 2)
Tumor_df$TumorLogFC <- TumorLogFC

# Establish differentially abundant lipids and metabolites and save their names for visualization
Tumor_df$diffexpressed[TumorPval < 0.05 & TumorFC < 0.8] <- "Downregulated"
Tumor_df$diffexpressed[TumorPval < 0.05 & TumorFC > 1.2] <- "Upregulated"
Tumor_df$delabel <- names(TumorAverage)

# Stroma regions ####
Stroma_df <- data.frame(StromaPadj, StromaPval, StromaFC)

# Make a column for differential abudance
Stroma_df$diffexpressed <- "Not significant"

# Log2 transform the FC
StromaLogFC <- log(StromaFC, base = 2)
Stroma_df$StromaLogFC <- StromaLogFC

# Establish differentially abundant lipids and metabolites and save their names for visualization
Stroma_df$diffexpressed[StromaPval < 0.05 & StromaFC < 0.8] <- "Downregulated"
Stroma_df$diffexpressed[StromaPval < 0.05 & StromaFC > 1.2] <- "Upregulated"
Stroma_df$delabel <- names(TumorAverage)

#### Differential analysis - Tumor vs Stroma ####
# Use the same data, except, now also exclude the sample with < 2% stroma from the tumor regions
TumorAverage <- TumorAverage[Inflam$stroma_perc > 2,]

# Initialize for the loop
Pval <- NA
FC <- NA

# Perform the t-tests and calculate P values and the fold changes (FC) for each lipid and metabolite
for(i in 1:295){
  # Use a paired t-test because we're comparing regions within the same samples
  Pval[i] <- t.test(TumorAverage[, i],StromaAverage[, i], paired = T)$p.val
  FC[i] <- mean(TumorAverage[, i]) / mean(StromaAverage[, i])
}

# Calculate the False Discovery Rate (FDR)
Padj <- p.adjust(Pval, method = "fdr")

# Create the data frame for the visualization of the results
DF <- data.frame(FC, Pval, Padj)

# Make a column for differential abudance
DF$diffexpressed <- "Not significant"

# Log2 transform the FC
DF$LogFC<-log(FC, base=2)

# Establish differentially abundant lipids and metabolites and save their names for visualization
DF$diffexpressed[Padj < 0.05 & DF$LogFC < -0.5] <- "Downregulated"
DF$diffexpressed[Padj < 0.05 & DF$LogFC > 0.5] <- "Upregulated"
DF$delabel<-names(TumorAverage)

#### Save the data ####
# Inflamed vs Non-inflamed regions
save(Tumor_df, Stroma_df, file = "./analysis_files/MALDI_Inflamed_vs_NonInflamed_regions.RData")

# Tumor vs Stroma regions
save(DF, file = "./analysis_files/MALDI_Tumor_vs_Stroma_regions.RData")
