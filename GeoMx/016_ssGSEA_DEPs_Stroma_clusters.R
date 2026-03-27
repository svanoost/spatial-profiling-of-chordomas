## Script for differential pathway enrichment analysis comparing inflamed with non-inflamed stroma segments
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(GeomxTools)
library(Biobase)
library(dplyr)
library(parallel)
library(GSVA)
library(GSEABase)
library(ComplexHeatmap)
library(circlize)

# Set working directory
master.location <- setwd(master.location)

# Read dataset
Chord_dataset = readRDS('./analysis_files/5_Chord_stroma_clusters.RData')

# Load WikiPathways list
own.WP <- getGmt("./input_files/Biologically_relevant_updated.gmt", sep = "\t")
names.WP <- names(own.WP)
for(i in 1:length(names.WP)){
  names.WP[i] <- strsplit(names.WP[i], split = "%")[[1]][1]
}
names(own.WP) <- names(names.WP)

#### Inflamed vs Non-inflamed stroma ####
# Perform the differential pathway enrichment analysis based on the segments
set.seed(1)
test_dataset = Chord_dataset
pData(test_dataset)$TEST = pData(test_dataset)$cluster

# Perform single sample Gene Set Enrichment Analysis (ssGSEA) per AOI
ssgsea_par <- ssgseaParam(exprData = assayDataElement(test_dataset,
                                                      elt = "log2_q3"),
                          geneSets = own.WP,
                          minSize = 5,
                          maxSize = 500)

ssgsea_results <- gsva(ssgsea_par)

# Make a Nanostring GeoMx Set from the ssGSEA results
geneSetObj.wp <-
  NanoStringGeoMxSet(assayData = ssgsea_results,
                     phenoData = AnnotatedDataFrame(pData(test_dataset)),
                     protocolData = protocolData(test_dataset),
                     featureType = "GeneSet",
                     check = FALSE)

# Use a linear mixed model (LMM) to correct for segments from the same sample
mixedOutmc <- mixedModelDE(geneSetObj.wp,
                           elt = "exprs",
                           modelFormula = ~ TEST + (1 | Tags),
                           groupVar = "TEST",
                           nCores = detectCores(),
                           multiCore = FALSE)

# Format results as data.frame
r_ssgsea_test <- do.call(rbind, mixedOutmc["lsmeans", ])
ssgsea_tests <- rownames(r_ssgsea_test)

r_ssgsea_test <- as.data.frame(r_ssgsea_test)
r_ssgsea_test$Contrast <- ssgsea_tests

r_ssgsea_test$Pathway <-
  unlist(lapply(colnames(mixedOutmc),
                rep, nrow(mixedOutmc["lsmeans", ][[1]])))

# Set the pathway names correctly
for(i in 1:nrow(r_ssgsea_test)){
  r_ssgsea_test[i, "Pathway"] <- strsplit(r_ssgsea_test[i, "Pathway"], split = "%")[[1]][1]
}

for(i in 1:nrow(ssgsea_results)){
  row.names(ssgsea_results)[i] <- strsplit(row.names(ssgsea_results)[i], split = "%")[[1]][1]
}

r_ssgsea_test$FDR <- p.adjust(r_ssgsea_test$`Pr(>|t|)`, method = "fdr")
r_ssgsea_test <- r_ssgsea_test[, c("Pathway", "Contrast", "Estimate",
                                   "Pr(>|t|)", "FDR")]

# Identify the top 10% of differentially enriched pathways (DEPs) for visualization
lmm_ssgsea_results <- r_ssgsea_test
es_thresh <- quantile(abs(lmm_ssgsea_results$Estimate), 0.9)

# Set Significance based on FDR < 0.05 and the percentile threshold
lmm_ssgsea_results$result_category = 'NS or NE'

lmm_ssgsea_results$result_category[
  lmm_ssgsea_results$FDR <= 0.05 &
    lmm_ssgsea_results$Estimate > es_thresh
] = 'Non-inflamed'

lmm_ssgsea_results$result_category[
  lmm_ssgsea_results$FDR <= 0.05 &
    lmm_ssgsea_results$Estimate < -es_thresh
] = 'Inflamed'

#### Generate the pathway enrichment heatmap ####
# Set the colors for annotation
segment_colors = c('#663300', '#f17a7b')
names(segment_colors) = c('TME', 'Tumor')

indiv_colors = c("#A349A4", "#FFFF33", "#E7298A", "#091833", "#1B9E77", "#D95F02", "#7570B3",  "#66A61E", "#8DD3C7",
                 "#9F000F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
names(indiv_colors) = unique(pData(Chord_dataset)$Tags)

cluster_colors <- c("#DC0018", "#2D69A9")
names(cluster_colors) <- c("2", "1") # Inflamed and Non-inflamed

nuclei_colors <- colorRamp2(c(0, 1500), c("white", "#5aff45"))

# Heatmap annotation for the columns
col_ann = HeatmapAnnotation(Segment = pData(Chord_dataset)$Segment,
                            Sample = pData(Chord_dataset)$Tags,
                            Cluster = pData(Chord_dataset)$cluster,
                            Nuclei = pData(Chord_dataset)$Nuclei,
                            col = list(Sample = indiv_colors, Nuclei = nuclei_colors,
                                    Segment = segment_colors, Cluster = cluster_colors),
                            annotation_legend_param = list(
                                    Sample = list(nrow = 2),
                                    Nuclei = list(direction = "horizontal")
                                    ))

# Make an object with the DEPs for the heatmap
DEPs = lmm_ssgsea_results[lmm_ssgsea_results$result_category == "Inflamed" | lmm_ssgsea_results$result_category == "Non-inflamed", "Pathway"]
names(DEPs) = lmm_ssgsea_results[lmm_ssgsea_results$result_category == "Inflamed" | lmm_ssgsea_results$result_category == "Non-inflamed", "result_category"]

# Calculate the Z-score
zscores_ssGSEA = ssgsea_results
for(i in 1:dim(zscores_ssGSEA)[1]){
  zscores_ssGSEA[i,] = (zscores_ssGSEA[i,] - mean(zscores_ssGSEA[i,])) / sd(zscores_ssGSEA[i,])
}

#### Save the data ####
save(DEPs, zscores_ssGSEA, col_ann, file = "./analysis_files/Fig_S6A_Heatmap.RData")
