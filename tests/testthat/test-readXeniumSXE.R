dir <- system.file(
  file.path("extdata", "Xenium_small"), 
  package="SpatialExperimentIO")
# dir <- "~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small"
dir <- "~/Downloads/BC_data/Xenium_rep2"
dir <- "~/Downloads/CRC_data/Xenium"

test_that("example data folders uniquely contains needed files", {
  expect_true("cell_feature_matrix.h5" %in% list.files(dir))
  # expect_true("cells.csv.gz" %in% list.files(dir)) 
  expect_true("cells.parquet" %in% list.files(dir)) 
  expect_true("cell_boundaries.parquet" %in% list.files(dir))
  expect_true("nucleus_boundaries.parquet" %in% list.files(dir))
  expect_true("transcripts.parquet" %in% list.files(dir)) 
  expect_true("experiment.xenium" %in% list.files(dir)) 
  
  expect_length(list.files(dir, "cell_feature_matrix.h5"), 1)
  expect_length(list.files(dir, "cells.parquet"), 1)
  expect_length(list.files(dir, "cell_boundaries.parquet"), 1)
  expect_length(list.files(dir, "nucleus_boundaries.parquet"), 1)
  expect_length(list.files(dir, "transcripts.parquet"), 1)
  expect_length(list.files(dir, "experiment.xenium"), 1)
})

test_that("needed files contains spatial columns of interest", {
  metadata <- read.csv(gzfile(file.path(dir, "cells.csv.gz")), header = TRUE)
  
  expect_true(all(c("x_centroid", "y_centroid") %in% colnames(metadata))) 
  expect_true(is.numeric(metadata$x_centroid))
  expect_true(is.numeric(metadata$y_centroid))
})

test_that("data are read correctly to SpatialExperiment class", {
  x <- readXeniumSXE(dir, 
                     returnType = "SPE",
                     countMatPattern = "cell_feature_matrix.h5",
                     metaDataPattern = "cells.parquet", # or cells.csv.gz 
                     coordNames = c("x_centroid", "y_centroid"), 
                     addExperimentXenium = TRUE,
                     altExps = NULL,
                     addParquetPaths = FALSE)
  
  expect_s4_class(x, "SpatialExperiment")
  expect_true(all(colnames(SpatialExperiment::spatialCoords(x)) == c("x_centroid", "y_centroid")))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(4, 6)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
  expect_true("experiment.xenium" %in% names(metadata(x)))
})

test_that("data are read correctly to SingleCellExperiment class", {
  x <- readXeniumSXE(dir, 
                     returnType = "SCE",
                     countMatPattern = "cell_feature_matrix.h5",
                     metaDataPattern = "cells.csv.gz", # back compatibility 
                     coordNames = c("x_centroid", "y_centroid"), 
                     addExperimentXenium = FALSE,
                     altExps = NULL,
                     addParquetPaths = TRUE,
                     loadTx = FALSE)
  
  expect_s4_class(x, "SingleCellExperiment")
  expect_true(all(c("x_centroid", "y_centroid") %in% colnames(SingleCellExperiment::colData(x))))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(4, 6)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
  expect_false("transcripts" %in% names(metadata(x)))
  expect_true("cell_boundaries" %in% names(metadata(x)))
  expect_true("nucleus_boundaries" %in% names(metadata(x)))
})

## Test altExp
# test_that("data are read correctly to SingleCellExperiment class", {
#   x <- readXeniumSXE(dir, 
#                      returnType = "SCE",
#                      countMatPattern = "cell_feature_matrix.h5",
#                      metaDataPattern = "cells.csv.gz", # back compatibility 
#                      coordNames = c("x_centroid", "y_centroid"), 
#                      addExperimentXenium = FALSE,
#                      altExps = c(negprobe="^NegControlProbe", antisense = "^antisense", negcode="^NegControlCodeword", blank = "^BLANK"),
#                      addParquetPaths = TRUE,
#                      loadTx = FALSE)
#   
#   expect_s4_class(x, "SingleCellExperiment")
#   expect_true(all(c("negprobe", "antisense", "negcode", "blank") %in% altExpNames(x)))
# })
