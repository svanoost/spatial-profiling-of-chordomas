## Script for loading the raw GeoMx data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(GeomxTools)

# Set working directory
master.location <- setwd(master.location)

# Get files' paths
ddc_files = dir(file.path(OG_data_dir, "dccs"), pattern=".dcc$", full.names=TRUE, recursive=TRUE)
pkc_files = dir(file.path(OG_data_dir, "pkcs"), pattern=".pkc$", full.names=TRUE, recursive=TRUE)
metadata_file = dir(file.path(OG_data_dir, "annotation"), pattern = "^[^~]", full.names=TRUE, recursive=TRUE)

# Read data into a NanoStringGeoMxSet object
Chord_dataset = readNanoStringGeoMxSet(dccFiles=ddc_files, pkcFiles=pkc_files,
                                       phenoDataFile=metadata_file, phenoDataSheet="Sheet1", phenoDataDccColName="Sample_ID",
                                       protocolDataColNames=c("aoi", "roi"), experimentDataColNames=c("panel"))

# Change column names of metadata
colnames(Chord_dataset@phenoData@data) = c('SlideName', 'ScanName', 'Segment', 'Area', 'Tags', 'SegmentDisplayName', 'Nuclei', 'ANN1',
                                        'ROICoordinateX', 'ROICoordinateY')

View(Chord_dataset@phenoData@data)

# Remove duplicate/ not important annotations:
Chord_dataset@phenoData@data = Chord_dataset@phenoData@data[-c(2,6,8)] # scanname, segmentdisplayname and ann1

#### Save the data ####
saveRDS(Chord_dataset, file = './input_files/1_Chord_dataset_raw.RData')
