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

# -------------
# Downsize to 10 genes and 9 cells for transcripts 
# -------------
path_to_tx <- "~/Downloads/Lung9_Rep1"
tx_cssv <- read.csv(file.path(path_to_tx, "Lung9_Rep1_tx_file.csv"))

# counts_test <- read.csv(file.path(cosmx_lung_p9s1_demo_path, "lung_p9s1_exprMat_file.csv"), row.names = FALSE)
ten_genes <- c("AATK","ABL1","ABL2","ACE","ACE2","ACKR1","ACKR3","ACKR4")
nine_cells <- 1:9
fov <- 10

tx_cssv_test <- tx_cssv %>%
  filter(fov == 10 & cell_ID %in% 1:9 & target %in% ten_genes)
  
write.csv(tx_cssv_test, file.path(cosmx_lung_p9s1_demo_path, "lung_p9s1_tx_file.csv"), row.names = FALSE) # only 1 row

# add to gitignore to mimic raw file download, but use locally for testing
arrow::write_parquet(tx_cssv_test, file.path(cosmx_lung_p9s1_demo_path, "tx_file.parquet"))  
# tx_cssv_test <- arrow::read_parquet(file.path(cosmx_lung_p9s1_demo_path, "tx_file.parquet")) # sanity


# -------------
# Downsize to 10 genes and 9 cells for fov position
# -------------
path_to_tx <- "~/Downloads/Lung9_Rep1"
fov_cssv <- read.csv(file.path(path_to_tx, "Lung9_Rep1_fov_positions_file.csv"))

fov_cssv_test <- fov_cssv %>%
  filter(fov == 10)

write.csv(fov_cssv_test, file.path(cosmx_lung_p9s1_demo_path, "lung_p9s1_fov_positions_file.csv"), row.names = FALSE) # only 1 row














