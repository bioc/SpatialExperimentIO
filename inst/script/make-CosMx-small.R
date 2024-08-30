################################################################################
# Script to create a small count matrix and a small coordinate files 
# from the raw download of the CosMx human lung cancer
# patient 9 slice 1.
# Yixing Dong, updated Aug 2024
################################################################################

# references:
# CosMx `_exprMat_file.csv` and `_metadata_file.csv` file were downloaded from
# \url{https://nanostring.com/resources/smi-ffpe-dataset-lung9-rep1-data/}

# in this script we subset the large raw necessary files into small ones.


# -------------
# Download data
# -------------

# Put the downloaded unzipped file into a folder. Make sure that two mandatory 
# files do exist. 

cosmx_lung_p9s1_path <- here::here("raw_data/cosmx_patient9slice1")
cosmx_folder <- list.files(cosmx_lung_p9s1_path, pattern = ".csv")
cosmx_folder


# -------------
# Read in raw data
# -------------

counts <- read.csv(file.path(cosmx_lung_p9s1_path, "Lung5_Rep1_exprMat_file.csv")) # 91992   982
meta <- read.csv(file.path(cosmx_lung_p9s1_path, "Lung5_Rep1_metadata_file.csv"))  # 91972    20

# -------------
# Downsize to 10 genes and 9 cells 
# -------------

library(dplyr)
counts_test <- counts %>%
  dplyr::filter(fov %in% c(10)) %>%
  slice(1:10) %>%
  dplyr::filter(cell_ID != 0) # take 9 cells

counts_test <- counts_test[, 1:10] # take 10 genes

meta_test <- meta %>%
  right_join(counts_test[, c("fov", "cell_ID")], by = c("fov", "cell_ID")) %>%
  select("fov", "cell_ID", "Area", "CenterX_local_px", "CenterY_local_px", "CenterX_global_px", "CenterY_global_px")

# -------------
# Save the small data 
# -------------

cosmx_lung_p9s1_demo_path <- "~/Desktop/SpatialExperimentIO/inst/extdata/CosMx_small"
write.csv(counts_test, file.path(cosmx_lung_p9s1_demo_path, "lung_p9s1_exprMat_file.csv"), row.names = FALSE)
write.csv(meta_test, file.path(cosmx_lung_p9s1_demo_path, "lung_p9s1_metadata_file.csv"), row.names = FALSE)





