################################################################################
# Script to create a small count matrix and a small coordinate files 
# by subsetting the raw seqFISH data provided in the closed Github issue #3 by 
# Github user methornton.
# Yixing Dong, updated Aug 2024
################################################################################

# As the original data was confidential, and only the data structure or column
# name is needed, I manually subsetted to the first 9 genes and 14 cells to 
# create the following data.

# Another 20GB example data from the official company website of Spatial Genomics
# can be found here \url{https://spatialgenomics.com/data/#kidney-data}


# -------------
# Make data
# -------------

cellxgene <- rbind(
  c("","ITGB1","NR4A1","THBS1","HIF1A","ID2","CD37","DUSP2","ETS1","STAT3"),
  c(1,0,0,0,0,0,0,0,0,0),
  c(2,0,0,0,0,0,0,0,0,0),
  c(3,0,0,0,0,0,0,0,0,0),
  c(4,0,0,0,1,0,0,0,0,0),
  c(5,0,0,0,0,0,0,0,0,0),
  c(6,0,0,0,0,0,0,0,0,0),
  c(7,1,0,0,0,0,0,0,0,0),
  c(8,0,0,0,0,0,0,0,0,0),
  c(9,0,0,0,0,0,0,0,0,0),
  c(10,0,0,0,0,0,0,0,0,0),
  c(11,0,0,0,0,0,0,0,0,0),
  c(12,0,0,0,0,0,0,0,0,0),
  c(13,0,0,0,0,0,0,0,0,0),
  c(14,0,0,0,0,0,0,0,0,0)
)

cellcoord <- rbind(
  c("label","area","center_x","center_y"),
  c(1,4293,23720,1512),
  c(2,6549,24498,1530),
  c(3,2670,23152,1514),
  c(4,4639,23537,1515),
  c(5,2284,23411,1517),
  c(6,1673,23659,1517),
  c(7,4943,23635,1544),
  c(8,4515,23457,1552),
  c(9,5645,23378,1562),
  c(10,5438,23212,1554),
  c(11,4677,23570,1569),
  c(12,1804,23512,1560),
  c(13,2171,23851,1551),
  c(14,3007,23045,1556)
)


# -------------
# Save the small data 
# -------------

seqfish_demo_path <- "~/Desktop/SpatialExperimentIO/inst/extdata/seqFISH_small"
write.csv(cellxgene, file.path(seqfish_demo_path, "TestSetSmall_CellxGene.csv"), row.names = FALSE)
write.csv(cellcoord, file.path(seqfish_demo_path, "TestSetSmall_CellCoordinates.csv"), row.names = FALSE)


# -------------
# Sanity check
# -------------

# install necessary packages

devtools::install_github("estellad/SpatialExperimentIO")

# load packages

library(SpatialExperimentIO)

# read the small data as SPE object

spe <- readSeqfishSXE(path, return_type = "SPE"); spe


