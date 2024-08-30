################################################################################
# Script to create a small count matrix and a small coordinate files 
# from the raw download of the STARmap PLUS mouse brain
# well 05.
# Yixing Dong, updated Aug 2024
################################################################################

# references:
# STARmap PLUS `raw_expression_pd.csv` and `spatial.csv` file were downloaded from
# \url{https://zenodo.org/records/8327576}

# in this script we subset the large raw necessary files into small ones.


# -------------
# Download data
# -------------
# Put the downloaded unzipped file into a folder. Make sure that two mandatory 
# files do exist. 
star_well05_path <- here::here("raw_data/star_well05")
star_folder <- list.files(star_well05_path, pattern = ".csv")
star_folder


# -------------
# Read in raw data
# -------------

counts <- read.csv(file.path(star_well05_path, "well05raw_expression_pd.csv")) # 91992   982
meta <- read.csv(file.path(star_well05_path, "well05_spatial.csv"))  # 91972    20


# -------------
# Downsize to 8 genes and 9 cells 
# -------------

library(dplyr)
meta_test <- meta %>%
  slice(1:10)  # take 9 cells 

cells <- meta_test$NAME[meta_test$NAME != "TYPE"]

counts_test <- counts %>%
  slice(1:8) # take 8 genes

counts_test <- counts_test[, colnames(counts_test) %in% c("GENE", cells)]


# -------------
# Save the small data 
# -------------

starmap_mousebrain_demo_path <- "~/Desktop/SpatialExperimentIO/inst/extdata/STARmapPLUS_small"
write.csv(counts_test, file.path(starmap_mousebrain_demo_path, "mockraw_expression_pd.csv"), row.names = FALSE)
write.csv(meta_test, file.path(starmap_mousebrain_demo_path, "mock_spatial.csv"), row.names = FALSE)


# -------------
# Sanity check
# -------------

# install necessary packages

devtools::install_github("estellad/SpatialExperimentIO")
remotes::install_github("lmweber/ggspavis")

# load packages

library(SpatialExperimentIO)
library(ggspavis)

# read the small data as SPE object, and plot it. 

spe <- readStarmapplusSXE(starmap_mousebrain_demo_path, return_type = "SPE")
plotSpots(spe, annotate = "Main_molecular_cell_type", in_tissue = NULL)



