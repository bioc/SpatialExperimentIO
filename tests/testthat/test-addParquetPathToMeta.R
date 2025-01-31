dir <- system.file(
  file.path("extdata", "CosMx_small"), 
  package="SpatialExperimentIO")

test_that("Able to add a file path for CosMx", {
  x <- readCosmxSXE(dirName = dir, addParquetPaths = FALSE)
  expect_true(length(metadata(x)) == 0)
  
  x <- addParquetPathToMeta(x, dirName = dir,
                            metaNames = "transcripts",
                            filePattern = "tx_file.csv")
  expect_true("transcripts" %in% names(metadata(x)))
})


dir <- system.file(
  file.path("extdata", "Xenium_small"), 
  package="SpatialExperimentIO")

test_that("Able to add a file path for Xenium", {
  x <- readXeniumSXE(dirName = dir, addParquetPaths = FALSE)
  expect_true(length(metadata(x)) == 1)
  
  x <- addParquetPathToMeta(x, dirName = dir, 
                            metaNames = "transcripts",
                            filePattern = "transcripts.parquet")
  expect_true("transcripts" %in% names(metadata(x)))
  
  x <- addParquetPathToMeta(x, dirName = dir, 
                            metaNames = "cell_boundaries",
                            filePattern = "cell_boundaries.parquet")
  expect_true("cell_boundaries" %in% names(metadata(x)))
  
  x <- addParquetPathToMeta(x, dirName = dir, 
                            metaNames = "nucleus_boundaries",
                            filePattern = "nucleus_boundaries.parquet")
  expect_true("nucleus_boundaries" %in% names(metadata(x)))
})
