
#' Example \code{SingleCellExperiment}
#' 
#' Contains pseudobulk data (mean expression per cell-type) from 11 different datasets. 
#' Mean expression matrices have been downsampled to 1,000/21,000 genes. 
#' 
#' 
#' @examples
#' \dontrun{
#' set.seed(2021)
#' pseudo_sce <- scKirby::ingest_data("/Users/schilder/Desktop/model_celltype_conservation/raw_data/scRNAseq/merged/merged.11datasets.h5ad", save_output = F)
#' SingleCellExperiment::reducedDim(pseudo_sce,"UMAP") <- data.frame(SummarizedExperiment::colData(pseudo_sce)[,c("UMAP.1","UMAP.2")])
#' pseudo_sce <- pseudo_sce[sample(1:nrow(pseudo_sce),1000), ] 
#' usethis::use_data(pseudo_sce, overwrite = T)
#' }
"pseudo_sce"




#' Example \code{Seurat}
#' 
#' Contains pseudobulk data (mean expression per cell-type) from 11 different datasets. 
#' Mean expression matrices have been downsampled to 1,000/21,000 genes. 
#' 
#' @examples
#' \dontrun{
#' set.seed(2021)
#' pseudo_seurat <- scKirby::ingest_data("/Users/schilder/Desktop/model_celltype_conservation/raw_data/scRNAseq/merged/merged.11datasets.h5ad", output_type="Seurat", save_output = F)
#' pseudo_seurat <- pseudo_seurat[sample(1:nrow(pseudo_seurat),1000), ]
#' usethis::use_data(pseudo_seurat, overwrite = T)
#' }
"pseudo_seurat"

