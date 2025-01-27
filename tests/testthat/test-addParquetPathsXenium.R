dir <- system.file(
  file.path("extdata", "Xenium_small"), 
  package="SpatialExperimentIO")
# dir <- "~/Desktop/SpatialExperimentIO/inst/extdata/Xenium_small"
# dir <- "~/Downloads/BC_data/Xenium_rep2"
# dir <- "~/Downloads/CRC_data/Xenium"

test_that("Able to add files independently for Xenium", {
  x <- readXeniumSXE(dirName = dir, addParquetPaths = FALSE)
  expect_true(length(metadata(x)) == 0)
  
  x <- addParquetPathsXenium(x, dirName = dir)
  expect_true("transcripts" %in% names(metadata(x)))
  expect_true("cell_boundaries" %in% names(metadata(x)))
  expect_true("nucleus_boundaries" %in% names(metadata(x)))
})
