#' @rdname readCosmxSXE
#' 
#' @title Load data from a Nanostring CosMx experiment
#' 
#' @description
#' Creates a \code{\link{SpatialExperiment}} from the downloaded unzipped CosMx  
#' directory for Nanostring CosMx spatial gene expression data.
#'
#' @param dirName a directory path to CosMx download that contains files of interest.
#' @param returnType option of \code{"SPE"} or \code{"SCE"}, stands for 
#' \code{SpatialExperiment} or \code{SingleCellExperiment} object. Default value \code{"SPE"}
#' @param countMatPattern a filename pattern for the count matrix. Default value is 
#' \code{"exprMat_file.csv"}, and there is no need to change.
#' @param metaDataPattern a filename pattern for the metadata .csv file that 
#' contains spatial coords. Default value is \code{"metadata_file.csv"}, and 
#' there is no need to change.
#' @param coordNames a vector of two strings specify the spatial coord names. 
#' Default value is \code{c("CenterX_global_px", "CenterY_global_px")}, and 
#' there is no need to change.
#' 
#' @details
#' The constructor assumes the downloaded unzipped CosMx folder has the following
#' structure, with two mandatory files:
#' CosMx_unzipped/optional_default_folder/ \cr
#' · | — *_exprMat_file.csv \cr
#' · | — *_metadata_file.csv \cr
#' 
#'
#' @return  a \code{\link{SpatialExperiment}} or a \code{\link{SingleCellExperiment}} object 
#' @export
#' 
#' @author Yixing Estella Dong
#'
#' @examples
#' # A relatively small data download can be from:
#' # https://nanostring.com/resources/smi-ffpe-dataset-lung9-rep1-data/
#' 
#' 
#' # A mock counts and mock metadata with spatial location generated for a 8 genes by 
#' # 9 cells object is in /extdata: 
#' 
#' cospath <- system.file(
#'   file.path("extdata", "CosMx_small"),
#'   package = "SpatialExperimentIO")
#'   
#' list.files(cospath)
#' 
#' # One of the following depending on your output (`SPE` or `SCE`) requirement.
#' cos_spe <- readCosmxSXE(dirName = cospath)
#' cos_sce <- readCosmxSXE(dirName = cospath, returnType = "SCE")
#' 
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment rowData counts colData
#' @importFrom methods as
#' @importFrom utils read.csv
readCosmxSXE <- function(dirName = dirName, 
                         returnType = "SPE",
                         countMatPattern = "exprMat_file.csv", 
                         metaDataPattern = "metadata_file.csv", 
                         coordNames = c("CenterX_global_px", "CenterY_global_px")){
  
  returnType <- match.arg(returnType, choices = c("SPE", "SCE"))
  
  ## Metadata sanity check 
  if(!any(file.exists(file.path(dirName, list.files(dirName, metaDataPattern))))){
    stop("CosMx metadata file does not exist in the directory. Expect 'metadata_file.csv' in `dirName`")
  }
  
  metadata_file <- file.path(dirName, list.files(dirName, metaDataPattern))
  if(length(metadata_file) > 1){
    stop("More than one metadata file possible with the provided pattern `metaDataPattern`")
  }
  
  ## Count matrix sanity check
  if(!any(file.exists(file.path(dirName, list.files(dirName, countMatPattern))))){
    stop("CosMx count matrix does not exist in the directory. Expect 'exprMat_file.csv' in `dirName`")
  }
  
  countmat_file <- file.path(dirName, list.files(dirName, countMatPattern))
  if(length(countmat_file) > 1){
    stop("More than one count matrix file possible with the provided pattern `countMatPattern`")
  }

  # Read in 
  countmat <- read.csv(countmat_file)
  metadata <- read.csv(metadata_file)
  
  # Count matrix   
  counts_ <- merge(countmat, metadata[, c("fov", "cell_ID")])
  counts_ <- subset(counts_, select = -c(fov, cell_ID))
  counts_ <- t(counts_)
  
  # rowData (does not exist)
  
  # metadata
  metadata <- merge(metadata, countmat[, c("fov", "cell_ID")])
  
  if(!all(coordNames %in% colnames(metadata))){
    stop("`coordNames` not in columns of `metaDataPattern`. For CosMx, expect c('CenterX_global_px', 'CenterY_global_px') in the columns of the metadata 'metadata_file.csv'. " )
  }
  
  colnames(counts_) <- rownames(metadata) <- seq_len(ncol(counts_))
  
  if(returnType == "SPE"){
  sxe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = as(counts_, "dgCMatrix")),
    # rowData = rowData,
    colData = metadata,
    spatialCoordsNames = coordNames)
  }else if(returnType == "SCE"){
    # construct 'SingleCellExperiment'
    sxe <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = as(counts_, "dgCMatrix")),
      colData = metadata
    )
  }
  
  return(sxe)
}
