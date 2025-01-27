dir <- system.file(
  file.path("extdata", "CosMx_small"), 
  package="SpatialExperimentIO")

# dir <- "~/Desktop/SpatialExperimentIO/inst/extdata/CosMx_small"
# dir <- "~/Downloads/Pancreas_CosMx_WTX"
# dir <- "~/Downloads/HumanFC"
# dir <- "~/Downloads/Lung9_Rep1"

test_that("Able to add files independently for CosMx", {
  x <- readCosmxSXE(dirName = dir, addParquetPaths = FALSE)
  expect_true(length(metadata(x)) == 0)
  
  x <- addParquetPathsCosMx(x, dirName = dir, loadPolygon = FALSE)
  expect_true("transcripts" %in% names(metadata(x)))
})
