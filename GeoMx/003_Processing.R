## Script for performing the QC for the GeoMx data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(GeomxTools)
library(NanoStringNCTools)
library(Biobase)
library(ggplot2)
library(uwot)
library(scales)

# Set working directory
master.location <- setwd(master.location)

setwd("I:/BEEN/Siddh/LCCO/Git_repositories/GEOMX_MALDI/GeoMx_repo")

# Read QC raw dataset:
Chord_dataset = readRDS('./analysis_files/2_Chord_dataset_QC.RData')

# Get modules used:
pkcs = annotation(Chord_dataset)
modules = gsub(".pkc", "", pkcs)

# Set the colors for annotations
segment_colors = c('#663300', '#f17a7b')
names(segment_colors) = c('TME', 'Tumor')

indiv_colors = c("#A349A4", "#FFFF33", "#E7298A", "#091833", "#1B9E77", "#D95F02", "#7570B3",  "#66A61E", "#8DD3C7",
                 "#9F000F", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
names(indiv_colors) = unique(pData(Chord_dataset)$Tags)

#### Create gene level count data ####
Chord_dataset = aggregateCounts(Chord_dataset)
dim(exprs(Chord_dataset)) #18 677 genes for 126 segments

#### Get the negative probes ####
negativeProbefData = subset(fData(Chord_dataset), CodeClass == "Negative")
neg_probes = unique(negativeProbefData$TargetName)

#### Check which segments to remove ####
# Calculate the limit of quantification (LOQ) per segment
cutoff = 2
minLOQ = 2
LOQ = data.frame(row.names = colnames(Chord_dataset))
for(module in modules) {
  vars = paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(Chord_dataset)))) {
    LOQ[, module] = pmax(minLOQ,
                         pData(Chord_dataset)[, vars[1]] * pData(Chord_dataset)[, vars[2]] ^ cutoff)
  }
}
pData(Chord_dataset)$LOQ <- LOQ

# Determine the number of genes detected in each segment based on the LOQ
LOQ_Mat = c()
for(module in modules) {
  ind = fData(Chord_dataset)$Module == module
  Mat_i = t(esApply(Chord_dataset[ind, ], MARGIN = 1,
                             FUN = function(x) {
                               x > LOQ[, module]
                               }))
  LOQ_Mat = rbind(LOQ_Mat, Mat_i)
}
LOQ_Mat = LOQ_Mat[fData(Chord_dataset)$TargetName, ]

pData(Chord_dataset)$GenesDetected = colSums(LOQ_Mat, na.rm=TRUE)
pData(Chord_dataset)$GeneDetectionRate = pData(Chord_dataset)$GenesDetected / nrow(Chord_dataset)
pData(Chord_dataset)$DetectionThreshold = cut(pData(Chord_dataset)$GeneDetectionRate,
                                                    breaks=c(0, 0.01, 0.05, 0.1, 0.15, 1),
                                                    labels=c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
View(pData(Chord_dataset))

# Check visually number of segments for each interval of gene detection rate
# Based on segments
ggplot(pData(Chord_dataset), aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = segment_colors) +
  labs(x = "Gene Detection Rate", y = "Number Segments", fill = "Segment")
# Based on samples
ggplot(pData(Chord_dataset), aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Tags)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = indiv_colors) +
  labs(x = "Gene Detection Rate", y = "Number Segments", fill = "Tags")

# Q3 of each segment should be higher than the negative probe count:
q3_stat_before = data.frame(row.names = colnames(exprs(Chord_dataset)),
                            Sample = colnames(exprs(Chord_dataset)),
                            QC_Flag = pData(Chord_dataset)[, 'QC_Flags'],
                            Q3 = unlist(apply(exprs(Chord_dataset), 2, quantile, 0.75, na.rm=TRUE)),
                            NegProbe = exprs(Chord_dataset)[neg_probes, ])
q3_stat_before$Q3_Flag = q3_stat_before$Q3 < q3_stat_before$NegProbe

ggplot(q3_stat_before, aes(x = NegProbe, y = Q3, color = Q3_Flag)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point(size = 2) + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  scale_color_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean", y = "Q3 Value")
# No Q3 is below the negative probe count

# Do a quick UMAP of all segments and see if flagged segments cluster apart:
# Update Flag info
pData(Chord_dataset)$DG1_Flags = pData(Chord_dataset)$DetectionThreshold == '<1%'
pData(Chord_dataset)$DG5_Flags = pData(Chord_dataset)$DetectionThreshold %in% c('<1%', '1-5%')
pData(Chord_dataset)$Q3_Flags = q3_stat_before[rownames(pData(Chord_dataset)), 'Q3_Flag']

# Normalize dataset
Chord_NORM_dataset = normalize(Chord_dataset, norm_method = "quant", desiredQuantile = .75, toElt = "q3_norm")

# Calculate UMAP with all genes
set.seed(1234)
umap_res = umap(t(log2(assayDataElement(Chord_NORM_dataset , elt = "q3_norm"))))

# Visualize UMAP
umap_df = data.frame(UMAP1 = umap_res[,1], UMAP2 = umap_res[,2],
                     QC_Flags = pData(Chord_NORM_dataset)[rownames(umap_res), 'QC_Flags'],
                     DG1_Flags = pData(Chord_NORM_dataset)[rownames(umap_res), 'DG1_Flags'],
                     DG5_Flags = pData(Chord_NORM_dataset)[rownames(umap_res), 'DG5_Flags'],
                     Q3_Flags = pData(Chord_NORM_dataset)[rownames(umap_res), 'Q3_Flags'],
                     Segment = pData(Chord_NORM_dataset)[rownames(umap_res), 'Segment'],
                     Tags = pData(Chord_NORM_dataset)[rownames(umap_res), 'Tags'],
                     Nuclei = pData(Chord_NORM_dataset)[rownames(umap_res), 'Nuclei'],
                     SlideName = pData(Chord_NORM_dataset)[rownames(umap_res), 'SlideName'])
# Overall QC
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = QC_Flags)) +
  geom_point(size = 3) + scale_colour_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme_bw()
# Detection threshold <1%
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = DG1_Flags)) +
  geom_point(size = 3) + scale_colour_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme_bw()
# Detection threshold <5%
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = DG5_Flags)) +
  geom_point(size = 3) + scale_colour_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme_bw()
# Q3
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Q3_Flags)) +
  geom_point(size = 3) + scale_colour_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme_bw()
# Annotated for segment
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Segment)) +
  geom_point(size = 2) + scale_colour_manual(values = segment_colors) +
  theme_bw()
#Annotated for sample
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Tags)) +
  geom_point(size = 2) + scale_colour_manual(values = indiv_colors) +
  theme_bw()

#### Filter segments ####
# Remove segments flagged in Q3 (none)
segments_to_keep = rownames(q3_stat_before)[!q3_stat_before$Q3_Flag]
Chord_filtered_dataset = Chord_dataset[, segments_to_keep]
dim(Chord_filtered_dataset) # 18677 genes 126 AOIs
dim(Chord_dataset) # 18677 genes 126 AOIs

#### Check which genes to remove ####
# Determine the number of segments in which genes were detected:
LOQ_Mat = LOQ_Mat[, colnames(Chord_filtered_dataset)]

# Globally
fData(Chord_filtered_dataset)$DetectedSegmentsGlobal = rowSums(LOQ_Mat, na.rm = TRUE)
fData(Chord_filtered_dataset)$DetectionRateGlobal =
  fData(Chord_filtered_dataset)$DetectedSegmentsGlobal / nrow(pData(Chord_filtered_dataset))

# Check visually number of genes detected by percentage of segments
# All segments
#    1%  = ~ 1 segments (out of 126)
#    5%  = ~ 6 segments (out of 126)
#    10% = ~ 13 segments (out of 126)
#    20% = ~ 25 segments (out of 126)
plot_detectGlobal = data.frame(Freq = factor(c('>1', '>5', '>10', '>20', '>30', '>50'), 
                                             levels = c('>1', '>5', '>10', '>20', '>30', '>50')))
plot_detectGlobal$Number = unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                                   function(x) {sum(fData(Chord_filtered_dataset)$DetectionRateGlobal >= x)}))
plot_detectGlobal$Rate = plot_detectGlobal$Number / nrow(fData(Chord_filtered_dataset))
rownames(plot_detectGlobal) = plot_detectGlobal$Freq

ggplot(plot_detectGlobal, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
                     vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                                high = "dodgerblue3", midpoint = 0.65,
                                limits = c(0,1),
                                labels = percent) +
  theme_bw() +
  scale_y_continuous(labels = percent, limits = c(0,1),
                              expand = expansion(mult = c(0, 0))) +
  labs(x = "Minimum % of Segments (Global)", y = "Genes Detected")
# Take 5% as a cut-off

#### Filter genes and redo all checks ####
# Filter genes
Chord_filtered_dataset = Chord_filtered_dataset[fData(Chord_filtered_dataset)$DetectionRateGlobal >= 0.05 |
                                            fData(Chord_filtered_dataset)$TargetName %in% neg_probes, ]
dim(Chord_filtered_dataset) # 8 401 genes and 126 AOIs

# Calculate the limit of quantification (LOQ) per segment
cutoff = 2
minLOQ = 2
LOQ = data.frame(row.names = colnames(Chord_filtered_dataset))
for(module in modules) {
  vars = paste0(c("NegGeoMean_", "NegGeoSD_"),
                module)
  if(all(vars[1:2] %in% colnames(pData(Chord_filtered_dataset)))) {
    LOQ[, module] = pmax(minLOQ,
                         pData(Chord_filtered_dataset)[, vars[1]] * pData(Chord_filtered_dataset)[, vars[2]] ^ cutoff)
  }
}
pData(Chord_filtered_dataset)$LOQ <- LOQ

# Determine the number of genes detected in each segment based on the LOQ
LOQ_Mat = c()
for(module in modules) {
  ind = fData(Chord_filtered_dataset)$Module == module
  Mat_i = t(esApply(Chord_filtered_dataset[ind, ], MARGIN = 1,
                             FUN = function(x) {
                               x > LOQ[, module]
                             }))
  LOQ_Mat = rbind(LOQ_Mat, Mat_i)
}
LOQ_Mat = LOQ_Mat[fData(Chord_filtered_dataset)$TargetName, ]

pData(Chord_filtered_dataset)$GenesDetected = colSums(LOQ_Mat, na.rm = TRUE)
pData(Chord_filtered_dataset)$GeneDetectionRate = pData(Chord_filtered_dataset)$GenesDetected / nrow(Chord_filtered_dataset)
pData(Chord_filtered_dataset)$DetectionThreshold = cut(pData(Chord_filtered_dataset)$GeneDetectionRate,
                                                             breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
                                                             labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))
View(pData(Chord_filtered_dataset))

# Check visually number of segments for each interval of gene detection rate
# Based on segments
ggplot(pData(Chord_filtered_dataset), aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = segment_colors) +
  labs(x = "Gene Detection Rate", y = "Number Segments", fill = "Segment")
# Based on samples
ggplot(pData(Chord_filtered_dataset), aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Tags)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = indiv_colors) +
  labs(x = "Gene Detection Rate", y = "Number Segments", fill = "Tags")

# Q3 of each segment should be higher than the negative probe count:
q3_stat_after = data.frame(row.names = colnames(exprs(Chord_filtered_dataset)),
                            Sample = colnames(exprs(Chord_filtered_dataset)),
                            QC_Flag = pData(Chord_filtered_dataset)[, 'QC_Flags'],
                            Q3 = unlist(apply(exprs(Chord_filtered_dataset), 2, quantile, 0.75, na.rm = TRUE)),
                            NegProbe = exprs(Chord_filtered_dataset)[neg_probes, ])
q3_stat_after$Q3_Flag = q3_stat_after$Q3 < q3_stat_after$NegProbe

ggplot(q3_stat_after, aes(x = NegProbe, y = Q3, color = Q3_Flag)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point(size = 2) + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  scale_color_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean", y = "Q3 Value")
# No segment with Q3 below negative probe count

# Do a quick UMAP of all segments and see if flagged segments cluster apart:
# Update Flag info
pData(Chord_filtered_dataset)$DG1_Flags = pData(Chord_filtered_dataset)$DetectionThreshold == '<1%'
pData(Chord_filtered_dataset)$DG5_Flags = pData(Chord_filtered_dataset)$DetectionThreshold %in% c('<1%', '1-5%')
pData(Chord_filtered_dataset)$Q3_Flags = q3_stat_after[rownames(pData(Chord_filtered_dataset)), 'Q3_Flag']

# Normalize dataset
Chord_filtered_NORM_dataset = normalize(Chord_filtered_dataset, norm_method = "quant", desiredQuantile = .75, toElt = "q3_norm")

# Calculate UMAP with all genes
umap_filt_res = uwot::umap(t(log2(assayDataElement(Chord_filtered_NORM_dataset , elt = "q3_norm"))))

# Visualize UMAP:
umap_df = data.frame(UMAP1 = umap_filt_res[,1], UMAP2 = umap_filt_res[,2],
                     QC_Flags = pData(Chord_filtered_NORM_dataset)[rownames(umap_filt_res), 'QC_Flags'],
                     DG1_Flags = pData(Chord_filtered_NORM_dataset)[rownames(umap_filt_res), 'DG1_Flags'],
                     DG5_Flags = pData(Chord_filtered_NORM_dataset)[rownames(umap_filt_res), 'DG5_Flags'],
                     Segment = pData(Chord_filtered_NORM_dataset)[rownames(umap_filt_res), 'Segment'],
                     Tags = pData(Chord_filtered_NORM_dataset)[rownames(umap_filt_res), 'Tags'],
                     SlideName = pData(Chord_filtered_NORM_dataset)[rownames(umap_filt_res), 'SlideName'])
# Overall QC
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = QC_Flags)) +
  geom_point(size = 3) + scale_colour_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme_bw()
# Detection threshold <1%
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = DG1_Flags)) +
  geom_point(size = 3) + scale_colour_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme_bw()
# Detection threshold <5%
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = DG5_Flags)) +
  geom_point(size = 3) + scale_colour_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'grey')) +
  theme_bw()
# Annotated for segment
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Segment)) +
  geom_point(size = 2) + scale_colour_manual(values = segment_colors) +
  theme_bw()
# Annotated for segment and sample
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Tags, shape = Segment)) +
  geom_point(size = 2) + scale_colour_manual(values = indiv_colors) + 
  scale_shape_manual(values = c('TME' = 1, 'Tumor' = 16))+
  theme_bw()
# Annotated for slide ID - batch effects?
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = SlideName)) +
  geom_point(size = 2) + 
  theme_bw()

# Log2 normalize the data
assayDataElement(object = Chord_filtered_NORM_dataset, elt = "log2_q3") = assayDataApply(Chord_filtered_NORM_dataset, 2, 
                                                                                     FUN = log, base = 2, elt = "q3_norm")
#### Save the data ####
# Filtered data
saveRDS(Chord_filtered_dataset, file = './analysis_files/3_Chord_dataset_FILTERED.RData')

# Filtered and normalized data
saveRDS(Chord_filtered_NORM_dataset, file = './analysis_files/4_Chord_dataset_NORM.RData')
