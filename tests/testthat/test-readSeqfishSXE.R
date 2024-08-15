dir <- system.file(
  file.path("extdata", "seqFISH_small"), 
  package="SpatialExperimentIO")

test_that("example data folders uniquely contains needed files", {
  expect_true("TestSetSmall_CellCoordinates.csv" %in% list.files(dir))
  expect_true("TestSetSmall_CellxGene.csv" %in% list.files(dir)) 
  
  expect_length(list.files(dir, "CellCoordinates.csv"), 1)
  expect_length(list.files(dir, "CellxGene.csv"), 1)
})

test_that("needed files contains spatial columns of interest", {
  metadata <- read.csv(file.path(dir, "TestSetSmall_CellCoordinates.csv"))
  
  expect_true(all(c("center_x", "center_y") %in% colnames(metadata))) 
  expect_true(is.numeric(metadata$center_x))
  expect_true(is.numeric(metadata$center_y))
})

test_that("data are read correctly to SpatialExperiment class", {
  x <- readSeqfishSXE(dirname = dir, 
                       return_type = "SPE",
                       countmatfpattern = "CellxGene.csv", 
                       metadatafpattern = "CellCoordinates.csv", 
                       coord_names = c("center_x", "center_y"))
  
  expect_s4_class(x, "SpatialExperiment")
  expect_true(all(colnames(SpatialExperiment::spatialCoords(x)) == c("center_x", "center_y")))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(9, 14)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
})

test_that("data are read correctly to SingleCellExperiment class", {
  x <- readSeqfishSXE(dirname = dir, 
                       return_type = "SCE",
                       countmatfpattern = "CellxGene.csv", 
                       metadatafpattern = "CellCoordinates.csv", 
                       coord_names = c("center_x", "center_y"))
  
  expect_s4_class(x, "SingleCellExperiment")
  expect_true(all(c("center_x", "center_y") %in% colnames(SingleCellExperiment::colData(x))))
  expect_false(is.null(rownames(x)))
  expect_false(is.null(colnames(x)))
  expect_true(all(dim(x) == c(9, 14)))
  expect_s4_class(SingleCellExperiment::counts(x), "dgCMatrix")
})
