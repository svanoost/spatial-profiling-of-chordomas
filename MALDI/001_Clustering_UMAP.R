## Script for generating the UMAP and clusters of the MALDI-MSI data
# R version 4.4.1

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(Seurat)
library(dplyr)
library(data.table)

# Set working directory
master.location <- setwd(master.location)

# Read the data (peak matrices with the Total Ion Counts (TIC) for all samples)
load("peakMat_merged.RData")

#### Prepare the data for UMAP and clustering ####
# Keep only columns with lipid & metabolite data
MyData <- peakMat_merge[, 2:296]

# load the identified lipid and metabolite names
peaks <- read.delim("./input_files/MALDI_peak_names.txt", header = T)$peaks

# Set the lipid and metabolite names as column names
colnames(MyData) <- peaks

# Select only the lipids and filter out all other molecules
lipids <- peaks[164:295]
# Filter out all unidentified molecules (starting with "X")
lipids <- lipids[!grepl("X", lipids)]
# Filter out molecules with "glu" in their name (UDP-glucose and oxidized glutathione)
lipids <- lipids[!grepl("glu", lipids)]
# Keep only lipids in the data for clustering and UMAP
MyData <- MyData[, colnames(MyData) %in% lipids]

# Set the Spot IDs as row names
rownames(MyData) = peakMat_merge$Spot_ID

# Normalize the TIC of the lipids
MyData <- MyData * 100
MyData_norm <- round(MyData, digits = 3)
MyData_norm <- as.data.frame(t(MyData_norm))

#### Perform UMAP and clustering using Seurat ####
# Create a Seurat object
Sobject <- CreateSeuratObject(counts = MyData_norm, project = "All samples combined")
# Log transform the data
Sobject <- NormalizeData(Sobject, normalization.method = "LogNormalize", scale.factor = 10)

# Check the variable features in the dataset and visualize
Sobject <- FindVariableFeatures(Sobject, selection.method = "vst", nfeatures = 41)
top10 <- head(VariableFeatures(Sobject), 10)
plot1 <- VariableFeaturePlot(Sobject, cols = c("black", "red"))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Perform the PCA
all.lipids <- rownames(x = Sobject)
Sobject <- ScaleData(object = Sobject, features = all.lipids)
Sobject <- RunPCA(object = Sobject, npcs = 20)

# Check how many PCs to use for clustering
ElbowPlot(object = Sobject, ndims = 20)

# Perform the clustering and UMAP
# Based on the Elbow plot you can change the dimension for the FindNeighbors
# Here we went with 11 PCs, resolution of 0.17, and 10 nearest neighbors after checking multiple variations
Sobject <-FindNeighbors(Sobject, dims = 1:11, reduction = "pca")
Sobject <- FindClusters(Sobject, resolution = 0.17, method = "igraph")
Sobject <- RunUMAP(object = Sobject, dims = 1:11, min.dist = 0.02, n.neighbors = 10)

# Extract the UMAP and cluster information from the Seurat object
UMAPimage <- data.frame(peakMat_merge[, 1:300], Sobject$seurat_clusters, Sobject@reductions$umap@cell.embeddings)
# Set the new column names
colnames(UMAPimage)[301:303] <- c("clusters", "UMAP1", "UMAP2")

#### Save the data ####
# Save the Seurat object
save(Sobject, file = "MALDI_Seurat_object_clusters.RData")

# Save the UMAP and clustering metadata using fwrite (makes for much faster saving and reading of the data)
fwrite(UMAPimage[, c("X", "L_ID", "x", "y", "xy", "clusters", "UMAP1", "UMAP2")],
       file = "MALDI_UMAP_clustering_metadata.txt", row.names = F, sep = "\t")
