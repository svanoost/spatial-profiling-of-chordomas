## Script for differential pathway enrichment analysis comparing tumor with stroma segments
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

# Set working directory
master.location <- setwd(master.location)

# Read dataset
Chord_dataset = readRDS('./analysis_files/4_Chord_dataset_NORM.RData')

# Load WikiPathways list
own.WP <- getGmt("./input_files/Biologically_relevant_WikiPathways.gmt", sep = "\t")
names.WP <- names(own.WP)
for(i in 1:length(names.WP)){
  names.WP[i] <- strsplit(names.WP[i], split = "%")[[1]][1]
}
names(own.WP) <- names(names.WP)

#### Tumor vs Stroma ####
# Perform the differential pathway enrichment analysis based on the segments
set.seed(1)
test_dataset = Chord_dataset
pData(test_dataset)$TEST = pData(test_dataset)$Segment

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

# Identify the top 10 differentially enriched pathways (DEPs) per segment
lmm_ssgsea_results <- r_ssgsea_test
lmm_ssgsea_results$result_category = 'NS or NE'

lmm_ssgsea_results <- lmm_ssgsea_results %>% arrange(Estimate)
lmm_ssgsea_results[1:10, "result_category"] <- "Tumor"
lmm_ssgsea_results <- lmm_ssgsea_results %>% arrange(-Estimate)
lmm_ssgsea_results[1:10, "result_category"] <- "Stroma"

# Create a separate data frame for the top DEPs for visualization
top_paths <- lmm_ssgsea_results[lmm_ssgsea_results$result_category != "NS or NE",]
top_paths$pathway_factor <- factor(top_paths$Pathway, levels = top_paths$Pathway)

#### Save the data ####
save(top_paths, file = "./analysis_files/Fig_1B_barplot.RData")
