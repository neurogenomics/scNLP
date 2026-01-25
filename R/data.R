

#' Example \code{Seurat}
#' 
#' Contains pseudobulk data (mean expression per cell-type) from 11 different datasets. 
#' Mean expression matrices have been downsampled to 1,000/21,000 genes. 
#' 
#' @examples
#' \dontrun{
#' set.seed(2021)
#' path <- "merged.11datasets.h5ad"
#' pseudo_seurat <- scKirby::ingest_data(path, output_type = "Seurat")
#' var_features <- FindVariableFeatures_split(
#'     pseudo_seurat, split.by = "batch",
#'     nfeatures = 4000, nfeatures_max = 2000
#' )
#' pseudo_seurat <- pseudo_seurat[var_features, ]
#' usethis::use_data(pseudo_seurat, overwrite = TRUE)
#' }
"pseudo_seurat"



#' Example \code{SingleCellExperiment}
#' 
#' Contains pseudobulk data (mean expression per cell-type) from 11 different datasets. 
#' Mean expression matrices have been downsampled to 1,000/21,000 genes. 
#' 
#' 
#' @examples
#' \dontrun{
#' set.seed(2021)
#' pseudo_sce <- scKirby::ingest_data(scNLP::pseudo_seurat)
#' umap_cols <- c("UMAP.1", "UMAP.2")
#' umap_df <- data.frame(
#'     SummarizedExperiment::colData(pseudo_sce)[, umap_cols]
#' )
#' SingleCellExperiment::reducedDim(pseudo_sce, "UMAP") <- umap_df
#' usethis::use_data(pseudo_sce, overwrite = TRUE)
#' }
"pseudo_sce"



