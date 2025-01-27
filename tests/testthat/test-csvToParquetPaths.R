dir <- system.file(file.path("extdata", "CosMx_small"),
                   package = "SpatialExperimentIO")
tx_csv_path <- file.path(dir, "lung_p9s1_tx_file.csv")

test_that("Able to convert .csv to .parquet for CosMx", {
  tx_parquet_path <- csvToParquetPaths(dir, filepath = tx_csv_path)
  
  expect_true(file.exists(tx_parquet_path))
})
