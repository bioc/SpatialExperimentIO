# ---------------------------
# Create metadata spreadsheet
# ---------------------------

# metadata for all datasets

df_all <- data.frame(
  Genome = NA, 
  SourceType = "FASTQ", 
  SourceVersion = NA, 
  Coordinate_1_based = NA, 
  DataProvider = NA, 
  Maintainer = "Yixing E. Dong <estelladong729@gmail.com>", 
  stringsAsFactors = FALSE
)


# metadata for individual datasets
df_Xenium_small <- cbind(
  df_all, 
  Title = "Xenium_small", 
  Description = paste0(
    "10x Genomics Xenium mock .h5 count matrix file and cells.csv.gz file ", 
    "downsized from the human breast cancer  replicate 1 in ",
    "Janesick et al. (2023)."), 
  SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7780153", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  BiocVersion = "3.19", 
  stringsAsFactors = FALSE
)


df_CosMx_small <- cbind(
  df_all, 
  Title = "CosMx_small", 
  Description =  paste0(
    "NanoString CosMx mock .csv count matrix file and .csv coordinates file ", 
    "downsized from the human NSCLC data (patient 9, slice 1). This dataset ", 
    "was previously released by NanoString on their website."), 
  SourceUrl = "https://nanostring.com/resources/smi-ffpe-dataset-lung9-rep1-data/", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  BiocVersion = "3.19", 
  stringsAsFactors = FALSE
)


df_MERSCOPE_small <- cbind(
  df_all, 
  Title = "MERSCOPE_small", 
  Description =  paste0(
    "Vizgen MERSCOPE mock .csv count matrix file and .csv coordinates file ", 
    "downsized from the human ovarian (patient 2, sample 1). This dataset was ", 
    "previously released by Vizgen on their website."), 
  SourceUrl = "https://console.cloud.google.com/storage/browser/vz-ffpe-showcase/HumanOvarianCancerPatient2Slice1;tab=objects?prefix=&forceOnObjectsSortingFiltering=false", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  BiocVersion = "3.19", 
  stringsAsFactors = FALSE
)


df_STARmapPLUS_small <- cbind(
  df_all, 
  Title = "STARmapPLUS_small", 
  Description =  paste0(
    "STARmap PLUS mock .csv count matrix file and .csv coordinates file ", 
    "downsized from the mouse brain (well 05), including annotations for ", 
    "cell type and tissue regions by Shi et al. (2023). "), 
  SourceUrl = "https://zenodo.org/records/8327576", 
  Species = "Mus musculus", 
  TaxonomyId = "10090", 
  BiocVersion = "3.19", 
  stringsAsFactors = FALSE
)


df_seqFISH_small <- cbind(
  df_all, 
  Title = "seqFISH_small", 
  Description = paste0(
    "Spatial Genomics seqFISH human unknown dataset obtained from Github", 
    " issue #3 by Github user methornton. "), 
  SourceUrl = "https://github.com/estellad/SpatialExperimentIO/issues/3/", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  BiocVersion = "3.19", 
  stringsAsFactors = FALSE
)


# combine and save as .csv spreadsheet file

df_combined <- rbind(
  df_Xenium_small, 
  df_CosMx_small, 
  df_MERSCOPE_small, 
  df_STARmapPLUS_small,
  df_seqFISH_small
)

write.csv(df_combined, file = "./inst/extdata/metadata.csv", row.names = FALSE)






