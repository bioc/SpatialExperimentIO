################################################################################
# Script to create a small count matrix and a small coordinate files 
# from the raw download of the Merscope human ovarian cancer
# patient 1 slice 2.
# Yixing Dong, updated Aug 2024
################################################################################

# references:
# Merscope `_cell_by_gene.csv` and `_cell_metadata.csv` file were downloaded from
# \url{https://console.cloud.google.com/storage/browser/vz-ffpe-showcase/HumanOvarianCancerPatient2Slice1;tab=objects?prefix=&forceOnObjectsSortingFiltering=false}

# in this script we subset the large raw necessary files into small ones.


# -------------
# Download data
# -------------

# Put the downloaded unzipped file into a folder. Make sure that two mandatory 
# files do exist. 

mer_lung_p1s2_path <- here::here("raw_data/mer_patient1slice2")
mer_folder <- list.files(mer_lung_p1s2_path, pattern = ".csv")
mer_folder


# -------------
# Read in raw data
# -------------

counts <- read.csv(file.path(mer_lung_p1s2_path, "OvarianP2S1_cell_by_gene.csv")) # 91992   982
meta <- read.csv(file.path(mer_lung_p1s2_path, "OvarianP2S1_cell_metadata.csv"))  # 91972    20


# -------------
# Downsize to 10 genes and 9 cells 
# -------------

library(dplyr)
meta_test <- meta[1:9, ] %>%
  select("X", "fov", "volume", "center_x", "center_y")

counts_test <- counts %>%
  mutate(X = cell) %>%
  right_join(meta_test[, c("X", "fov")], by = "X") %>%
  select(-c(X, fov)) 

counts_test <- counts_test[1:10, ]


# -------------
# Save the small data 
# -------------

mer_ovarian_p1s2_demo_path <- "~/Desktop/SpatialExperimentIO/inst/extdata/MERSCOPE_small"
write.csv(counts_test, file.path(mer_ovarian_p1s2_demo_path, "ovarian_p1s2_cell_by_gene.csv"), row.names = FALSE)
write.csv(meta_test, file.path(mer_ovarian_p1s2_demo_path, "ovarian_p1s2_cell_metadata.csv"), row.names = FALSE)


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

spe <- readMerscopeSXE(mer_ovarian_p1s2_demo_path, return_type = "SPE")
plotSpots(spe, annotate = "volume", in_tissue = NULL)

