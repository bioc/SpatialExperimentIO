#' Add CosMx-related parquet paths to metadata for transcripts, polygon, or cell/nucleus boundaries.
#' 
#' @param sxe a SPE or SCE Xenium object to add parquet to `metadata(sxe)`.
#' @param dirName the directory that stores the transcripts/polygon 
#' .csv or .parquet files.
#' @param loadTx to load path to transcripts parquet to \code{metadata(sxe)}or not. 
#' Default is TRUE. 
#' @param txMetaNames names to add to slots in \code{metadata(sxe)[["name"]]}. 
#' The number of `txMetaNames` should equal to number of file detected in `dirName` 
#' with `txPattern`. Can have multiple, such as \code{c("transcripts", "transcripts1")}. 
#' Default is \code{"transcripts"}.
#' @param txPattern .csv or .parquet (if you have previous converted) pattern 
#' of transcript file in `dirName`. Can have multiple, such as \code{c("tx_file.csv", "tx_file1.csv")}. 
#' Default value is \code{"tx_file.csv"}.
#' @param loadPolygon to load path to polygons parquet to \code{metadata(sxe)} or not. 
#' Default is TRUE.
#' @param polygonMetaNames names to add to slots in `metadata(sxe)$`. The number of 
#' `polygonMetaNames` should equal to number of file detected in `dirName` with `polygonPattern`.
#' Can have multiple. Can have multiple, such as \code{c("polygons", "polygons1")}. 
#' Default is \code{"transcripts"}.
#' @param polygonPattern .csv or .parquet (if you have previous converted) pattern
#'  of polygons file in the `dirName`. Can have multiple, such as \code{c("polygons.csv", "polygons1.csv")}. 
#' Default value is \code{"polygons.csv"}. 
#' 
#' @author Yixing Estella Dong
#' 
#' @return a SPE or SCE Xenium object with parquet paths added to metadata
#' @export 
#' 
#' @examples 
#' cospath <- system.file(file.path("extdata", "CosMx_small"), 
#'                        package = "SpatialExperimentIO")
#' 
#' sxe <- readCosmxSXE(dirName = cospath, addParquetPaths = FALSE)
#' sxe <- addParquetPathsCosMx(sxe, dirName = cospath, loadPolygon = FALSE)
#'
addParquetPathsCosMx <- function(sxe,
                                 dirName, 
                                 loadTx = TRUE,
                                 txMetaNames = "transcripts",
                                 txPattern = "tx_file.csv", 
                                 loadPolygon = TRUE,
                                 polygonMetaNames = "polygons",
                                 polygonPattern = "polygons.csv"){
  
  # Transcripts
  if(loadTx) sxe <- addParquetPathToMeta(sxe, dirName = dirName,
                                         metaNames = txMetaNames,
                                         filePattern = txPattern)
  
  # Polygon
  if(loadPolygon) sxe <- addParquetPathToMeta(sxe, dirName = dirName,
                                              metaNames = polygonMetaNames,
                                              filePattern = polygonPattern)
  
  return(sxe)
}


#' Add Xenium-related parquet paths to metadata for transcripts or cell/nucleus boundaries.
#' 
#' @param sxe a SPE or SCE Xenium object to add parquet to `metadata(sxe)`.
#' @param dirName the directory that stores the transcripts/cell_boundaries/nucleus_boundaries 
#' .parquet files.
#' @param loadTx to load path to transcripts parquet to \code{metadata(sxe)} or not. 
#' Default is FALSE. 
#' @param txMetaNames names to add to slots in \code{metadata(sxe)[["name"]]}. 
#' The number of `txMetaNames` should equal to number of file detected in `dirName` 
#' with `txPattern`. Can have multiple, such as \code{c("transcripts", "transcripts1")}. 
#' Default is \code{"transcripts"}.
#' @param txPattern .parquet pattern of transcript file in `dirName`. Can have multiple, 
#' such as \code{c("transcripts.parquet", "transcripts1.parquet")}. Default value is 
#' \code{"transcripts.parquet"}.
#' @param loadCellBound to load path to cell boundaries parquet to \code{metadata(sxe)} or not. 
#' Default is FALSE. 
#' @param cellBoundMetaNames names to add to slots in \code{metadata(sxe)[["name"]]}. 
#' The number of `cellBoundMetaNames` should equal to number of file detected in `dirName` 
#' with `cellBoundPattern`. Can have multiple, such as \code{c("cell_boundaries", "cell_boundaries1")}. 
#' Default is \code{"cell_boundaries"}.
#' @param cellBoundPattern .parquet pattern of cell boundaries file in `dirName`. Can have multiple, 
#' such as \code{c("cell_boundaries.parquet", "cell_boundaries1.parquet")}. Default value is 
#' \code{"cell_boundaries.parquet"}.
#' @param loadNucBound to load path to nucleus boundaries parquet to \code{metadata(sxe)} or not. 
#' Default is FALSE. 
#' @param NucBoundMetaNames names to add to slots in \code{metadata(sxe)[["name"]]}. 
#' The number of `NucBoundMetaNames` should equal to number of file detected in `dirName` 
#' with `NucBoundPattern`. Can have multiple, such as \code{c("nucleus_boundaries", "nucleus_boundaries1")}. 
#' Default is \code{"nucleus_boundaries"}.
#' @param NucBoundPattern .parquet pattern of nucleus boundaries file in `dirName`. Can have multiple, 
#' such as \code{c("nucleus_boundaries.parquet", "nucleus_boundaries1.parquet")}. Default value is 
#' \code{"nucleus_boundaries.parquet"}.
#' 
#' @author Yixing Estella Dong
#' 
#' @return a SPE or SCE Xenium object with parquet paths added to metadata
#' @export 
#'
#' @examples 
#' xepath <- system.file(file.path("extdata", "Xenium_small"),
#'                       package = "SpatialExperimentIO")
#' 
#' sxe <- readXeniumSXE(dirName = xepath, addParquetPaths = FALSE)
#' sxe <- addParquetPathsXenium(sxe, dirName = xepath)
#'
addParquetPathsXenium <- function(sxe, 
                                  dirName,
                                  loadTx = TRUE,
                                  txMetaNames = "transcripts",
                                  txPattern = "transcripts.parquet", 
                                  loadCellBound = TRUE,
                                  cellBoundMetaNames = "cell_boundaries",
                                  cellBoundPattern = "cell_boundaries.parquet", 
                                  loadNucBound = TRUE,
                                  NucBoundMetaNames = "nucleus_boundaries",
                                  NucBoundPattern = "nucleus_boundaries.parquet"){
  
  # Transcripts
  if(loadTx) sxe <- addParquetPathToMeta(sxe, dirName = dirName,
                                         metaNames = txMetaNames,
                                         filePattern = txPattern)
  
  # Cell Boundaries
  if(loadCellBound) sxe <- addParquetPathToMeta(sxe, dirName = dirName,
                                                metaNames = cellBoundMetaNames,
                                                filePattern = cellBoundPattern)
  
  # Nucleus Boundaries
  if(loadNucBound) sxe <- addParquetPathToMeta(sxe, dirName = dirName,
                                               metaNames = NucBoundMetaNames,
                                               filePattern = NucBoundPattern)
  
  return(sxe)
}





