## Script for performing the QC for the GeoMx data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(GeomxTools)
library(Biobase)

# Set working directory
master.location <- setwd(master.location)

# Read raw dataset:
Chord_dataset = readRDS('./input_files/1_Chord_dataset_raw.RData')

# Shift genes with 0 counts to 1
Chord_dataset = shiftCountsOne(Chord_dataset, useDALogic=TRUE)

#### Perform segment QC ####
# Check quality based on cutoffs and flag segments that do not pass
QC_params = list(minSegmentReads = 1000, # Minimum number of reads (1000)
                 percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
                 percentStitched = 80,   # Minimum % of reads stitched (80%)
                 percentAligned = 80,    # Minimum % of reads aligned (80%)
                 percentSaturation = 50, # Minimum sequencing saturation (50%)
                 minNegativeCount = 1,   # Minimum negative control counts (1)
                 maxNTCCount = 12000,    # Maximum counts observed in NTC well (12000)
                 minNuclei = 100,         # Minimum # of nuclei estimated (100)
                 minArea = 1250)         # Minimum segment area (1250)
Chord_dataset = setSegmentQCFlags(Chord_dataset, qcCutoffs=QC_params)
Chord_dataset@protocolData@data$QCFlags$LowNuclei = Chord_dataset@phenoData@data[rownames(Chord_dataset@protocolData@data), 
                                                                                 'Nuclei'] < 100
View(Biobase::protocolData(Chord_dataset)[["QCFlags"]])
# 7 segments have low alignment, 2 segments have low stitched

# To which individuals and segments do these non-passing segments belong to?
QCFlags = protocolData(Chord_dataset)[["QCFlags"]]
not_pass_segments = rownames(QCFlags)[rowSums(QCFlags) > 0]
Chord_dataset@phenoData@data$QC_Flags = FALSE
Chord_dataset@phenoData@data[not_pass_segments, 'QC_Flags'] = TRUE

View(Chord_dataset@phenoData@data)
# L3382 - TME; L3356 - TME; L3356 - Tumor; 2x L3770 - TME; L3519 - TME; L3519 - Tumor
# These will be excluded when they appear as outliers in the next step of the QC.
# If not, they will be kept for further analysis.

#### Perform probe QC ####
# Check probe quality based on cutoffs and flag probes that do not pass:
Chord_dataset = setBioProbeQCFlags(Chord_dataset,
                                   qcCutoffs=list(minProbeRatio=0.1, percentFailGrubbs = 20),
                                   removeLocalOutliers=TRUE)

# Check number of probes that passed or were flagged at global or local level:
ProbeQCResults = fData(Chord_dataset)[["QCFlags"]]
qc_df = data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                   Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                   Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                               & !ProbeQCResults$GlobalGrubbsOutlier))

# 3. Only 21 probes failed at local level. These will be excluded.
Chord_dataset = subset(Chord_dataset,
                       fData(Chord_dataset)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
                       fData(Chord_dataset)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

#### Save the data ####
saveRDS(Chord_dataset, file = './analysis_files/2_Chord_dataset_QC.RData')
