## Script for differential gene expression analysis comparing tumor with stroma segments
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(GeomxTools)
library(Biobase)
library(dplyr)
library(parallel)

# Set working directory
master.location <- setwd(master.location)

# Read dataset
Chord_dataset = readRDS('./analysis_files/4_Chord_dataset_NORM.RData')

#### Tumor vs Stroma ####
# Perform the differential gene expression analysis based on the segments
set.seed(1)
test_dataset = Chord_dataset
pData(test_dataset)$TEST = pData(test_dataset)$Segment

# Use the linear mixed model (LMM) that takes into account that there are multiple segments per sample
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
r_test$result_category[r_test$Estimate > 0.7] = 'Stroma'
r_test$result_category[r_test$Estimate < (-0.7)] = 'Tumor'
r_test$result_category[r_test$FDR > 0.05] = 'NS or NE'

# Identify the top 10 DEGs for each segment for visualization in volcano plot
r_test$top_genes <- "no"

r_test <- r_test %>% arrange(Estimate)
r_test[1:10, "top_genes"] <- "yes"
r_test <- r_test %>% arrange(-Estimate)
r_test[1:10, "top_genes"] <- "yes"

#### Generate the DEG heatmap ####
# Set the colors for annotations
segment_colors = c('#663300', '#f17a7b')
names(segment_colors) = c('TME', 'Tumor')

pheno_colors <- c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
names(pheno_colors) <- c("Stromal cells", "Plasma cells", "Stromal / Myeloid", "TBXT target", "Myeloid cells")

indiv_colors = c("#A349A4", "#FFFF33", "#E7298A", "#091833", "#1B9E77", "#D95F02", "#7570B3",  "#66A61E", "#8DD3C7",
                 "#9F000F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
names(indiv_colors) = unique(pData(Chord_dataset)$Tags)

# Heatmap annotation for the columns - AOIs
col_ann = HeatmapAnnotation(Segment = pData(Chord_dataset)$Segment,
                            Sample = pData(Chord_dataset)$Tags,
                            col = list(Sample = indiv_colors,
                                       Segment = segment_colors))

# Load annotation for the differentially expressed genes (DEGs)
annotation <- read.delim("./input_files/Tumor_Stroma_Clusters.txt")
annotation$Gene <- gsub("\\?", "-", annotation$Gene)
row.names(annotation) = annotation$Gene

# Make an object with the DEGs for the heatmap
DEGs = r_test[r_test$result_category == "Stroma" | r_test$result_category == "Tumor", "Gene"]
names(DEGs) = r_test[r_test$result_category == "Stroma" | r_test$result_category == "Tumor", "result_category"]

# Calculate the Z-score
zscores = assayDataElement(Chord_dataset, elt = 'log2_q3')
for(i in 1:dim(zscores)[1]){
  zscores[i,] = (zscores[i,] - mean(zscores[i,])) / sd(zscores[i,])
}

# Heatmap annotation for the rows - genes
annotation = annotation[row.names(zscores[DEGs,]),]
row_ann = rowAnnotation(Cluster = annotation$Cluster,
                        col = list(Cluster = pheno_colors))

#### Save the data ####
# DEG results Tumor vs Stroma
write.table(r_test, file = "./output_files/DEGs_Tumor_vs_Stroma_007_cutoff.txt", 
            sep = "\t", row.names = FALSE)

# Data for the heatmap - Supplementary Figure 3
save(col_ann, row_ann, DEGs, zscores, file = "./analysis_files/Suppl_Fig_3_Heatmap.RData")

# Data for the volcano plot - Figure 1A
save(r_test, file = "./analysis_files/Fig_1A_Volcanoplot.RData")
