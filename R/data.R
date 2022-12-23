

#' Example \code{Seurat}
#' 
#' Contains pseudobulk data (mean expression per cell-type) from 11 different datasets. 
#' Mean expression matrices have been downsampled to 1,000/21,000 genes. 
#' 
#' @examples
#' \dontrun{
#' set.seed(2021)
#' pseudo_seurat <- scKirby::ingest_data("/Volumes/Steelix/model_celltype_conservation/merged/merged.11datasets.h5ad", output_type="Seurat", save_output = F)
#' ### Find the most variable features within each dataset  
#' var_features <- FindVariableFeatures_split(pseudo_seurat, split.by = "batch", nfeatures = 4000, nfeatures_max = 2000)
#' pseudo_seurat <- pseudo_seurat[var_features, ]
#' usethis::use_data(pseudo_seurat, overwrite = T)
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
#' pseudo_sce <- scKirby::ingest_data(scNLP::pseudo_seurat, save_output = F)
#' SingleCellExperiment::reducedDim(pseudo_sce,"UMAP") <- data.frame(SummarizedExperiment::colData(pseudo_sce)[,c("UMAP.1","UMAP.2")])
#' usethis::use_data(pseudo_sce, overwrite = T)
#' }
"pseudo_sce"



