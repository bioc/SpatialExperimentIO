################################################################################
# Script to create a small .h5 count matrix and a small cells.csv.gz files 
# from the raw download of the Xenium 1 human breast cancer
# in 10x Xenium paper by Janesick et al. (2023)
# Yixing Dong, updated Aug 2024
################################################################################

# references:
# XeniumOutputBundle .zip file was downloaded from 
# \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7780153}

# in this script we subset the large raw necessary files into small ones.

# -------------
# Download data
# -------------

# In the 10X Xenium paper, the Xenium sample 1 (2 replicates) data was accompanied 
# by consecutive slices of Chromium and Visium replicates. 

# We first prepare the data folder required by reader function `SpatialExperimentIO::readXeniumSXE()`. 
# Please click to download [Xenium Sample 1 Replicate 1 outs](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7780153&format=file&file=GSM7780153%5FXenium%5FFFPE%5FHuman%5FBreast%5FCancer%5FRep1%5Fouts%2Ezip). 
# Here we placed the downloaded objects in `here::here("/raw_data/xenium_sample1rep1/")`. 

xenium_out_folder <- list.files(here::here("raw_data/xenium_sample1rep1"), pattern = ".zip")
xenium_out_folder


# We now unzip `spatial.tar.gz` to its same directory. 

unzip(file.path(here::here("raw_data/xenium_sample1rep1"), xenium_out_folder),
      exdir = here::here("raw_data/xenium_sample1rep1"))


# Have a look at all the files in the newly unzipped `/outs` folder. The only 
# files you need are `cell_feature_matrix.h5` or `cell_feature_matrix`, and 
# `cells.csv.gz`. 

list.files(here::here("raw_data/xenium_sample1rep1/outs"))

# Sanity check that the files you need are in the `/outs` folder.

all(c("cell_feature_matrix.h5", "cell_feature_matrix", "cells.csv.gz") %in% 
      list.files(here::here("raw_data/xenium_sample1rep1/outs")))


# -------------
# Read raw .h5 file
# -------------

# First, we need to install the `rhdf5` loader package that would 
# read the .h5 object: 

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("rhdf5")

library(rhdf5)

h5file <- "~/Desktop/SampleData/archive/Xenium_Preprint_Data/Xenium/outs/cell_feature_matrix.h5"
h5ls(file = h5file)


# -------------
# Downsize to 4 genes, 6 cells small .h5 for each sub-level of .h5
# -------------

h5_barcodes = h5read(h5file, "/matrix/barcodes")
h5_barcodes_new <- as(head(h5_barcodes), class(h5_barcodes)) #1 2 3 4 5 6 # take only 6 cells

h5_data = h5read(h5file, "/matrix/data")
# matrix(c(0, 2:6, 0, 8:10, 0, 0, 0, 14:15, 0, 17:18, 0, 20, 0, 2:3, 0), nrow = 4, ncol = 6) # something like this
# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices
h5_data_new <- c(2, 3, 4, 5, 6, 8, 9, 10, 14, 15, 17, 18, 20, 2, 3) # int 1 1 1 3 2 4  # Nonzero UMI counts in column-major order


## Features
h5_all_tag_keys = h5read(h5file, "/matrix/features/_all_tag_keys")
h5_all_tag_keys_new <- h5_all_tag_keys # "genome"

h5_feature_type = h5read(h5file, "/matrix/features/feature_type")
h5_feature_type_new <- head(h5_feature_type, 4) # "Gene Expression" # take only 4 genes

h5_genome = h5read(h5file, "/matrix/features/genome")
h5_genome_new <- head(h5_genome, 4) # "Unknown" # take only 4 genes

h5_id = h5read(h5file, "/matrix/features/id")
h5_id_new <- head(h5_id, 4) # "ENSG000121270" # take only 4 genes

h5_name = h5read(h5file, "/matrix/features/name")
h5_name_new <- head(h5_name, 4) # "ABCC11" # take only 4 genes

# matrix
h5_indices = h5read(h5file, "/matrix/indices")
h5_indices_new <- h5_data_new <- c(1, 2, 3, 0, 1, 3, 0, 1, 1, 2, 0, 1, 3, 1, 2) # int 33 36 67  # Zero-based row index of corresponding element in `data`

h5_indptr = h5read(h5file, "/matrix/indptr")
h5_indptr_new <- c(0, 3, 6, 8, 10, 13, 15) # int 33 36 67 # take 6 cells, # Zero-based index into data / indices of the start of each column, 
# i.e., the data corresponding to each barcode sequence
# https://stackoverflow.com/questions/52299420/scipy-csr-matrix-understand-indptr

h5_shape = h5read(h5file, "/matrix/shape")
h5_shape_new <- h5_shape # 541 167780
h5_shape_new[1] <- 4
h5_shape_new[2] <- 6


# -------------
# How orignal .h5 looks like
# -------------

#               group          name       otype  dclass      dim
# 0                 /        matrix   H5I_GROUP                 
# 1           /matrix      barcodes H5I_DATASET  STRING   167780
# 2           /matrix          data H5I_DATASET INTEGER 10640759
# 3           /matrix      features   H5I_GROUP                 
# 4  /matrix/features _all_tag_keys H5I_DATASET  STRING        1
# 5  /matrix/features  feature_type H5I_DATASET  STRING      541
# 6  /matrix/features        genome H5I_DATASET  STRING      541
# 7  /matrix/features            id H5I_DATASET  STRING      541
# 8  /matrix/features          name H5I_DATASET  STRING      541
# 9           /matrix       indices H5I_DATASET INTEGER 10640759
# 10          /matrix        indptr H5I_DATASET INTEGER   167781
# 11          /matrix         shape H5I_DATASET INTEGER        2


# -------------
# Format and write my small .h5
# -------------

smallh5 <- "~/Desktop/SampleData/SpatialExperimentIOData/Xenium/cell_feature_matrix.h5"

h5createFile(smallh5)
h5createGroup(smallh5,"matrix")
h5createGroup(smallh5,"matrix/features")
h5ls(smallh5)

h5write(h5_barcodes_new, smallh5,"matrix/barcodes")
h5write(h5_data_new, smallh5,"matrix/data")

h5write(h5_all_tag_keys_new, smallh5,"matrix/features/_all_tag_keys")
h5write(h5_feature_type_new, smallh5,"matrix/features/feature_type")
h5write(h5_genome_new, smallh5,"matrix/features/genome")
h5write(h5_id_new, smallh5,"matrix/features/id")
h5write(h5_name_new, smallh5,"matrix/features/name")

h5write(h5_indices_new, smallh5,"matrix/indices")
h5write(h5_indptr_new, smallh5,"matrix/indptr")
h5write(h5_shape_new, smallh5,"matrix/shape")
h5ls(smallh5)


# -------------
# Check if raw and small .h5 can be read as a SingleCellExperiment
# -------------

# Now try to read it as SCE

library(SingleCellExperiment)
# chromh5 <- "~/Desktop/SampleData/archive/Xenium_Preprint_Data/Chromium_FFPE_Human_Breast_Cancer_Chromium_FFPE_Human_Breast_Cancer/count_sample_filtered_feature_bc_matrix.h5"

sce <- DropletUtils::read10xCounts(h5file, type = "HDF5", col.names = TRUE)
# sce <- DropletUtils::read10xCounts(chromh5, type = "HDF5", col.names = TRUE)
scetest <- SingleCellExperiment::SingleCellExperiment(assays = assays(sce))
counts(scetest)[34:38,]

sce <- DropletUtils::read10xCounts(smallh5, type = "HDF5", col.names = TRUE)
scetest <- SingleCellExperiment::SingleCellExperiment(assays = assays(sce))

data_check <- counts(sce)


# -------------
# Create small metadata
# -------------

# install.packages("dplyr")
library(dplyr)

# Read un .gz ed cell.csv

cells <- read.csv("~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small/cells.csv")

# Subset to 6 cells

cells <- cells %>% slice(1:6)

# Save the subset

write.csv(cells, "~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small/cells.csv")

# Now manually zipped it to have `cells.csv.gz` in ./inst/extdata/Xenium_small/


# -------------
# Sanity check
# -------------

# install necessary packages

devtools::install_github("estellad/SpatialExperimentIO")
remotes::install_github("lmweber/ggspavis")

# load packages

library(SpatialExperimentIO)
library(ggspavis)

# read the small .h5 and cells.csv.gz as SPE object, and plot it. 

spe <- readXeniumSXE(dirName = "~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small/", returnType = "SPE")
plotSpots(spe, annotate = "total_counts", in_tissue = NULL)


sce <- DropletUtils::read10xCounts("~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small/cell_feature_matrix.h5", 
                                   type = "HDF5", col.names = TRUE)

rownames(sce) <- rowData(sce)$Symbol

four_genes <- c("ABCC11", "ACTA2", "ACTG2", "ADAM9")
six_cells <- 1:6


# -------------
# Create small nucleus boundaries
# -------------
filepaths <- "~/Downloads/BC_data/Xenium_rep1/"
nucboundname <- "nucleus_boundaries.parquet"
cellboundname <- "cell_boundaries.parquet"
txname <- "transcripts.parquet"
cellsname <- "cells.parquet"

xe_small_path <- "~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small"

nucbound <- arrow::read_parquet(file.path(filepaths, nucboundname))
nucbound_test <- nucbound %>% filter(cell_id %in% six_cells)
arrow::write_parquet(nucbound_test, file.path(xe_small_path, nucboundname))  

cellbound <- arrow::read_parquet(file.path(filepaths, cellboundname))
cellbound_test <- cellbound %>% filter(cell_id %in% six_cells)
arrow::write_parquet(cellbound_test, file.path(xe_small_path, cellboundname))  

tx <- arrow::read_parquet(file.path(filepaths, txname), as_data_frame=FALSE)
tx <- tx |>
  mutate(feature_name=as.character(feature_name)) |>
  mutate(transcript_id=as.numeric(transcript_id)) |>
  as.data.frame()
tx_test <- tx %>% filter(cell_id %in% six_cells & feature_name %in% four_genes)
arrow::write_parquet(tx_test, file.path(xe_small_path, txname))  

cells <- arrow::read_parquet(file.path(filepaths, cellsname))
cells_test <- cells %>% filter(cell_id %in% six_cells)
arrow::write_parquet(cells_test, file.path(xe_small_path, cellsname))  

## Sanity
# extdatapath <- "~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small/"
# tx <- arrow::read_parquet(file.path(extdatapath, txname))
# cells <- arrow::read_parquet(file.path(extdatapath, cellsname))
# cells <- arrow::read_parquet(file.path(extdatapath, cellsname))

