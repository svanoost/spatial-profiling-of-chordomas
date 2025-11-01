## Script for clustering of tumor segments only
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(NanoStringNCTools)
library(Biobase)
library(parallel)
library(dplyr)
library(uwot)
library(bluster)

# Set working directory
master.location <- setwd(master.location)

# Read dataset
Chord_dataset = readRDS('./analysis_files/4_Chord_dataset_NORM.RData')
Chord_tumor <- Chord_dataset[, pData(Chord_dataset)$Segment == 'Tumor']

# Coefficient of variation (CV) function
calc_CV = function(x){sd(x) / mean(x)}

#### UMAP of tumor segments ####
# Find genes with top coefficient of variation (CV)
# Calculate CV values
CV_data = assayDataApply(Chord_tumor, elt = "log2_q3", MARGIN = 1, calc_CV)

# Save CV values:
saveRDS(CV_data, './analysis_files/CV_data_tumor.RData')

# Get genes in the top 20% of the CV values:
GOI = names(CV_data)[CV_data > quantile(CV_data, 0.8)]

# UMAP of all segments:
set.seed(1)
umap_all = umap(t(assayDataElement(Chord_tumor[GOI, ], elt = "log2_q3")), n_neighbors = 6)

# Save UMAP
write.csv(umap_all, './output_files/1_UMAP_TumorSegments.csv')

# Reload UMAP
umap_all <- read.csv('./output_files/1_UMAP_TumorSegments.csv')
row.names(umap_all) <- umap_all$X

# Make a UMAP data frame
umap_df = data.frame(UMAP1 = umap_all[,2], UMAP2 = umap_all[,3],
                     Segment = pData(Chord_tumor)[rownames(umap_all), 'Segment'],
                     Individual = pData(Chord_tumor)[rownames(umap_all), 'Tags'],
                     Location = pData(Chord_tumor)[rownames(umap_all), 'Location'],
                     SlideName = pData(Chord_tumor)[rownames(umap_all), 'SlideName'],
                     Nuclei = pData(Chord_tumor)[rownames(umap_all), 'Nuclei'])

#### Graph-based clustering of tumor segments ####
set.seed(1)
clust_graphbased = clusterRows(t(assayDataElement(Chord_tumor[GOI, ], elt = "log2_q3")),
                                        SNNGraphParam(k = 6, cluster.fun = 'louvain', 
                                                      cluster.args = list(resolution = 0.4)))
umap_df$snn <- clust_graphbased
umap_df$snn <- as.character(umap_df$snn)

#### Differential gene expression analysis ####
# Add the graph-based clusters to the GeoMx dataset
pData(Chord_tumor)$cluster <- clust_graphbased

# Perform the differential analysis, similar to that of Tumor vs Stroma
set.seed(1)
test_dataset = Chord_tumor
pData(test_dataset)$TEST = pData(test_dataset)$cluster

# Use the LMM
mm_res = mixedModelDE(test_dataset,
                      elt = "log2_q3",
                      modelFormula = ~ TEST + (1 | Tags),
                      groupVar = "TEST",
                      nCores = detectCores(),
                      multiCore = FALSE)

r_test = do.call(rbind, mm_res["lsmeans", ])
tests = rownames(r_test)
r_test = as.data.frame(r_test)
r_test$Contrast = tests
r_test$Gene = unlist(lapply(colnames(mm_res), rep, nrow(mm_res["lsmeans", ][[1]])))
r_test$FDR = p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
r_test = r_test[, c("Gene", "Contrast", "Estimate", "Pr(>|t|)", "FDR")]

# Set the thresholds for significance
r_test$result_category = 'NS or NE'
r_test$result_category[r_test$Estimate > 1] = 'Non-inflamed'
r_test$result_category[r_test$Estimate < (-1)] = 'Inflamed'
r_test$result_category[r_test$FDR > 0.05] = 'NS or NE'

# Identify the top 20 DEGs for the inflamed cluster for the volcano plot
r_test <- r_test %>% arrange(Estimate)
r_test$top_genes <- "no"
r_test[1:20, "top_genes"] <- "yes"

#### Generate the DEG heatmap ####
# Set colors for annotation
segment_colors = c('#663300', '#f17a7b')
names(segment_colors) = c('TME', 'Tumor')

indiv_colors = c("#A349A4", "#FFFF33", "#E7298A", "#091833", "#1B9E77", "#D95F02", "#7570B3",  "#66A61E", "#8DD3C7",
                 "#9F000F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
names(indiv_colors) = unique(pData(Chord_tumor)$Tags)

cluster_colors <- c("#1B9E77", "#D95F02")
names(cluster_colors) <- c("2", "1")

# Heatmap annotation for the columns
col_ann = HeatmapAnnotation(Segment = pData(test_dataset)$Segment,
                            Sample = pData(test_dataset)$Tags,
                            Cluster = pData(test_dataset)$cluster,
                             col = list(Sample = indiv_colors, Segment = segment_colors, 
                                      Cluster = cluster_colors))

# Make an object with the DEGs for the heatmap
DEGs = r_test[r_test$result_category == "Inflamed" | r_test$result_category == "Non-inflamed", "Gene"]
names(DEGs) = r_test[r_test$result_category == "Inflamed" | r_test$result_category == "Non-inflamed", "result_category"]

# Calculate the Z-score
zscores = assayDataElement(test_dataset, elt = 'log2_q3')
for(i in 1:dim(zscores)[1]){
  zscores[i,] = (zscores[i,] - mean(zscores[i,])) / sd(zscores[i,])
}

#### Save the data ####
# DEG results for tumor clusters
write.table(r_test, file = "./output_files/DEGs_Tumor_1_cutoff_graph_based_clusters.txt", 
            sep = "\t", row.names = FALSE)

# GeoMx dataset with only tumor segments and their clusters
saveRDS(Chord_tumor, file = "./analysis_files/5_Chord_tumor_clusters.RData")

# Data for the heatmap - Supplementary Figure 6
save(col_ann, DEGs, zscores, file = "./analysis_files/Suppl_Fig_6_Heatmap.RData")

# Data for the UMAPs - Figure 1E
save(umap_df, file = "./analysis_files/Fig_1E_UMAPs.RData")

# Data for the volcano plot - Figure 1F
save(r_test, file = "./analysis_files/Fig_1F_Volcanoplot.RData")
