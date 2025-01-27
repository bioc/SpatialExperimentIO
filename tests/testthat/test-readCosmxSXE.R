dir <- system.file(
  file.path("extdata", "CosMx_small"), 
  package="SpatialExperimentIO")
# dir <- "~/Desktop/SpatialExperimentIO/inst/extdata/CosMx_small"
dir <- "~/Downloads/Pancreas_CosMx_WTX"
# dir <- "~/Downloads/HumanFC"
# dir <- "~/Downloads/Lung9_Rep1"

test_that("example data folders uniquely contains needed files", {
  expect_true("lung_p9s1_exprMat_file.csv" %in% list.files(dir))
  expect_true("lung_p9s1_metadata_file.csv" %in% list.files(dir)) 
  expect_true("lung_p9s1_fov_positions_file.csv" %in% list.files(dir))
  expect_true("lung_p9s1_tx_file.csv" %in% list.files(dir))
  
  # parquet
  expect_true("tx_file.parquet" %in% list.files(dir))
  
  expect_length(list.files(dir, "exprMat_file.csv"), 1)
  expect_length(list.files(dir, "metadata_file.csv"), 1)
  expect_length(list.files(dir, "fov_positions_file.csv"), 1)
  expect_length(list.files(dir, "tx_file.csv"), 1)
  
  # parquet
  expect_length(list.files(dir, "tx_file.parquet"), 1)
})

test_that("needed files contains spatial columns of interest", {
  metadata <- read.csv(file.path(dir, "lung_p9s1_metadata_file.csv"))
  
  expect_true(all(c("CenterX_global_px", "CenterY_global_px") %in% colnames(metadata))) 
  expect_true(is.numeric(metadata$CenterX_global_px))
  expect_true(is.numeric(metadata$CenterY_global_px))
})

test_that("needed fov position files contains columns of interest, with fov position in colData", {
  metadata <- read.csv(file.path(dir, "lung_p9s1_fov_positions_file.csv"))
  
  expect_true(all(c("fov", "x_global_px", "y_global_px") %in% colnames(metadata))) 
  expect_true(is.numeric(metadata$fov))
  expect_true(is.numeric(metadata$x_global_px))
  expect_true(is.numeric(metadata$y_global_px))
})

test_that("data are read correctly to SpatialExperiment class", {
  x <- readCosmxSXE(dirName = dir, 
                    returnType = "SPE",
                    countMatPattern = "exprMat_file.csv", 
                    metaDataPattern = "metadata_file.csv", 
                    coordNames = c("CenterX_global_px", "CenterY_global_px"),
                    loadFovPos = TRUE,
                    fovPosPattern = "fov_positions_file.csv",
                    altExps = NULL,
                    addParquetPaths = FALSE)
  
  expect_s4_class(x, "SpatialExperiment")
  expect_true(all(colnames(SpatialExperiment::spatialCoords(x)) == c("CenterX_global_px", "CenterY_global_px")))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(8, 9)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
  expect_true(all(c("fov", "x_global_px", "y_global_px") %in% colnames(x@colData))) 
})

test_that("data are read correctly to SingleCellExperiment class, no fov position in colData", {
  x <- readCosmxSXE(dirName = dir, 
                    returnType = "SCE",
                    countMatPattern = "exprMat_file.csv", 
                    metaDataPattern = "metadata_file.csv", 
                    coordNames = c("CenterX_global_px", "CenterY_global_px"),
                    loadFovPos = FALSE,
                    fovPosPattern = "fov_positions_file.csv",
                    altExps = NULL,
                    addParquetPaths = FALSE)
  
  expect_s4_class(x, "SingleCellExperiment")
  expect_true(all(c("CenterX_global_px", "CenterY_global_px") %in% colnames(SingleCellExperiment::colData(x))))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(8, 9)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
  expect_false(all(c("fov", "x_global_px", "y_global_px") %in% colnames(x@colData))) 
})

test_that("data can be read in, even no polygon file, had loadPolygon = FALSE, and transcripts read in to metadata(()", {
  x <- readCosmxSXE(dirName = dir, 
                    returnType = "SPE",
                    countMatPattern = "exprMat_file.csv", 
                    metaDataPattern = "metadata_file.csv", 
                    coordNames = c("CenterX_global_px", "CenterY_global_px"),
                    loadFovPos = TRUE,
                    fovPosPattern = "fov_positions_file.csv",
                    altExps = NULL,
                    addParquetPaths = TRUE,
                    loadPolygon = FALSE)
  expect_true(metadata(x)$transcripts == file.path(dir, "lung_p9s1_tx_file.parquet"))
  expect_true(is.null(metadata(x)$polygons))
})

test_that("data can be read in, directly read parquet of transcripts", {
  x <- readCosmxSXE(dirName = dir, 
                    returnType = "SPE",
                    countMatPattern = "exprMat_file.csv", 
                    metaDataPattern = "metadata_file.csv", 
                    coordNames = c("CenterX_global_px", "CenterY_global_px"),
                    loadFovPos = TRUE,
                    fovPosPattern = "fov_positions_file.csv",
                    altExps = NULL,
                    addParquetPaths = TRUE,
                    loadPolygon = FALSE,
                    txPattern = "tx_file.parquet")
  expect_true(metadata(x)$transcripts == file.path(dir, "lung_p9s1_tx_file.parquet"))
  # expect_true(metadata(x)$transcripts == file.path(dir, "tx_file.parquet"))
})

test_that("data can be read in, directly read parquet of transcripts multiple", {
  x <- readCosmxSXE(dirName = dir, 
                    returnType = "SPE",
                    countMatPattern = "exprMat_file.csv", 
                    metaDataPattern = "metadata_file.csv", 
                    coordNames = c("CenterX_global_px", "CenterY_global_px"),
                    loadFovPos = TRUE,
                    fovPosPattern = "fov_positions_file.csv",
                    altExps = NULL,
                    addParquetPaths = TRUE,
                    loadPolygon = FALSE,
                    txMetaNames = c("transcripts", "transcripts2"),
                    txPattern = c("tx_file.csv", "tx_file.parquet"))
  expect_true(metadata(x)$transcripts == file.path(dir, "lung_p9s1_tx_file.parquet"))
  expect_true(metadata(x)$transcripts2 == file.path(dir, "lung_p9s1_tx_file.parquet"))
})

## Test altExp
# test_that("data can be read in, directly read parquet of transcripts", {
#   x <- readCosmxSXE(dirName = dir, 
#                     returnType = "SPE",
#                     countMatPattern = "exprMat_file.csv", 
#                     metaDataPattern = "metadata_file.csv", 
#                     coordNames = c("CenterX_global_px", "CenterY_global_px"),
#                     loadFovPos = TRUE,
#                     fovPosPattern = "fov_positions_file.csv",
#                     altExps = c(negprobe="^Neg", falsecode="^Sys"),
#                     addParquetPaths = FALSE)
#   
#   expect_true(all(c("negprobe", "falsecode") %in% altExpNames(x)))
# })




