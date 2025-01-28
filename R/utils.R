#' Sanity check if one and only file with the specified name pattern exists in 
#' the data download directory, and return the file path to .csv
#' Used for count matrix and metadata only, as they require unique. 
#'
#' @param tech Name of technology. Defined at the beginning of the function. e.g. "CosMx"
#' @param filetype File type to do sanity check. e.g. "metadata"
#' @param expectfilename Expected file pattern name for this file type. e.g. "metadata_file.csv"
#' @param dirName Directory to the data download. 
#' @param filepatternvar The file pattern variable. e.g. "metaDataPattern"
#'
#' @author Yixing Estella Dong
#' 
#' @return a path to a unique file of count matrix or colData.
#'
#' @examples 
#' \dontrun{
#' dir <- system.file(file.path("extdata", "CosMx_small"),
#'                    package = "SpatialExperimentIO")
#' countmat_file <- SpatialExperimentIO:::.sanityCheck(tech = "CosMx", filetype = "count matrix",
#'                               expectfilename = "`exprMat_file.csv`",
#'                               dirName = dir,
#'                               filepatternvar = "exprMat_file.csv")
#'}
.sanityCheck <- function(tech, filetype, expectfilename, dirName, filepatternvar){
  if(!any(file.exists(file.path(dirName, list.files(dirName, filepatternvar))))){
    stop(paste(tech, filetype , "file does not exist in the directory. Expect", expectfilename, "in", "`dirName`"))
  }
  
  path_to_file <- file.path(dirName, list.files(dirName, filepatternvar))
  if(length(path_to_file) > 1){
    stop(paste("More than one", filetype, "file possible with the provided pattern"))
  }
  
  return(path_to_file)
}



#' If transcripts or polygon is expected to be loaded, write a parquet file to 
#' the current data download (if not already), and return the file path to .parquet
#'
#' @param dirName current directory of data download
#' @param filepath path to transcripts or polygons csv
#'
#' @author Yixing Estella Dong
#'
#' @return a path to .parquet
#' @export 
#'
#' @examples 
#' dir <- system.file(file.path("extdata", "CosMx_small"),
#'                    package = "SpatialExperimentIO")
#' tx_csv_path <- file.path(dir, "lung_p9s1_tx_file.csv")
#' tx_parquet_path <- csvToParquetPaths(dirName, filepath = tx_csv_path)
#' 
#' @importFrom data.table fread
#' @importFrom arrow write_parquet
#' 
csvToParquetPaths <- function(dirName, filepath = "tx_csv_path"){
  parquet_path <- paste0(gsub(".csv", "", filepath), ".parquet")
  if(!file.exists(parquet_path)) write_parquet(as.data.frame(fread(filepath)), parquet_path)

  return(parquet_path)
}



#' Add parquet paths to metadata for transcripts, polygon, or cell/nucleus boundaries.
#'
#' @param sxe a SPE or SCE object to add parquet to `metadata(sxe)`.
#' @param dirName the directory that stores the transcripts/polygon/cell_boundaries 
#' .csv or .parquet files.
#' @param metaNames  a vector of names to `metadata(sxe)[[]]`. The length must
#' match number of files detected with filePattern provided. 
#' \code{e.g. c("transcripts", "transcripts1.csv")}.
#' @param filePattern a vector of file patterns to search in the current directory. 
#' e.g. \code{c("tx_file.csv", "tx_file1.csv")}.
#'
#' @return a SPE or SCE object with parquet paths added to metadata
#' @export 
#'
#' @examples
#' dir <- system.file(file.path("extdata", "CosMx_small"),
#'                    package = "SpatialExperimentIO")
#' sxe <- readCosmxSXE(dir)
#' sxe <- addParquetPathToMeta(sxe,
#'                             dirName = dir,
#'                             metaNames = "transcripts",
#'                             filePattern = "tx_file.parquet")
#' 
#' @importFrom purrr walk2
#' @importFrom S4Vectors metadata metadata<-
#' 
addParquetPathToMeta <- function(sxe, 
                                 dirName = dirName,
                                 metaNames = "transcripts",
                                 filePattern = "tx_file.csv"){
  if(!all(grepl(".csv|.parquet", filePattern))){
    stop("Require each `*Pattern` have a '.csv' or '.parquet' file extension.")
  }
  
  fileswpat <- unlist(lapply(filePattern, list.files, path = dirName))
  
  if(length(metaNames) != length(fileswpat)){
    stop("Number of metadata slot names and number of detect files with `*Pattern` have different lengths")
  }
  
  if(length(fileswpat) == 0){
    stop("File with `*Pattern` does not exist in the directory.")
  }else{
    if(all(grepl(".csv", fileswpat))){           # all csv, no parquet
      tx_files <- file.path(dirName, fileswpat)
      parquet_paths <- lapply(tx_files, csvToParquetPaths, dirName = dirName)
    }else if(all(grepl(".parquet", fileswpat))){ # only parquet
      parquet_paths <- file.path(dirName, fileswpat[grepl(".parquet", fileswpat)])
    }else{                                       # a mix of parquet & csv
      tx_files <- file.path(dirName, fileswpat[grepl(".csv", fileswpat)])
      parquet_paths_csv <- lapply(tx_files, csvToParquetPaths, dirName = dirName)
      
      parquet_paths <- file.path(dirName, fileswpat[grepl(".parquet", fileswpat)])
      parquet_paths <- append(parquet_paths_csv, parquet_paths)
    }
  }
  
  parquet_paths <- unlist(parquet_paths)
  
  walk2(metaNames, parquet_paths, function(name, path) {
    metadata(sxe)[[name]] <<- path
  })
  
  return(sxe)
}



