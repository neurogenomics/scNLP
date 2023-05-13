#' Run standardized \pkg{Seurat} pipeline
#' 
#' Run \pkg{Seurat} pipeline on \pkg{Seurat} object 
#' or raw \code{counts} and \code{meta.data}.\cr
#' Automatically performs 
#' \describe{
#' \item{\code{FindVariableFeatures}}{Variable feature selection}
#' \item{\code{NormalizeData}}{Data normalization}
#' \item{\code{RunPCA}}{PCA}
#' \item{\code{RunUMAP}}{UMAP}
#' \item{\code{FindNeighbors}}{K-nearest neighbors}
#' \item{\code{FindClusters}}{Clustering}
#' }
#' @inheritDotParams scKirby::process_seurat
#' @returns A preprocessed \link[Seurat]{Seurat} object.
#' 
#' @export
#' @import Seurat
#' @import future
#' @examples 
#' X <- scNLP::pseudo_sce@assays@data$raw@seed[,seq_len(100)]
#' obj <- seurat_pipeline(obj = X)
seurat_pipeline <- function(...){
  .Deprecated("scKirby::process_seurat")  
  scKirby::process_seurat(...)
}



