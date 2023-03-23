
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
#' @param cluster_reduction Recompute neighbors graph based on UMAP 
#' to get clusters that best reflect UMAP space.
#' For this same reason, only cluster in two dimensions, 
#' because this is the view we most often use.
#' That said, this may reduce the generalisability of these clusters/graph.
#' @returns A preprocessed \link[Seurat]{Seurat} object.
#' 
#' @export
#' @import Seurat
#' @import future
#' @examples 
#' X <- scNLP::pseudo_sce@assays@data$raw@seed[,seq_len(100)]
#' obj <- seurat_pipeline(obj = X)
seurat_pipeline <- function(obj = NULL,
                            meta.data = NULL,
                            nfeatures = 2000,
                            vars.to.regress = NULL,
                            dims = seq_len(50),
                            add_specificity = FALSE,
                            assay_name="RNA",
                            default_assay = NULL,
                            n.components = 2L,
                            log_norm = FALSE,
                            cluster_reduction = "umap",
                            workers = 1,
                            max_mem = 8000*1024^2,
                            seed = 2020){
  requireNamespace("Seurat") 
  set.seed(seed)
  future::plan(strategy = "multicore", workers = workers)
  options(future.globals.maxSize = max_mem)
  
  if(!methods::is(obj,"Seurat")){
    seurat <- Seurat::CreateSeuratObject(counts = obj,  
                                         meta.data = meta.data,
                                         assay = assay_name)
  } else { 
    seurat <- obj
  } 
  #### Create new assay using EWCE-style specificity ####
  if(isTRUE(add_specificity)){
    print("+ Calculating specificity as a new assay...") 
    spec1 <- calc_specificity(as.matrix(Seurat::GetAssayData(seurat)))
    seurat[['specificity']] <- Seurat::CreateAssayObject(counts = spec1) 
  }
  
  ### Set assay ###
  if(!is.null(default_assay)) Seurat::DefaultAssay(seurat) <- default_assay 
  #### Select variable features ####
  if(is.null(nfeatures)) nfeatures <- nrow(seurat)
  seurat <- Seurat::FindVariableFeatures(seurat, 
                                         nfeatures = nfeatures)
  seurat <- Seurat::NormalizeData(seurat)
  if(isTRUE(log_norm)) logged <- Seurat::LogNormalize(seurat)
  seurat <- Seurat::ScaleData(seurat,  
                              vars.to.regress = vars.to.regress)
  
  
  seurat <- Seurat::RunPCA(seurat)
  seurat <- Seurat::FindNeighbors(seurat)
  seurat <- Seurat::RunUMAP(seurat,
                            dims = dims,
                            n.components = n.components, 
                            return.model = TRUE) 
  if(!isFALSE(cluster_reduction)){
    seurat <- Seurat::FindNeighbors(seurat,
                                    reduction = cluster_reduction, 
                                    dims = seq_len(n.components))
    seurat <- Seurat::FindClusters(seurat,
                                   reduction = cluster_reduction)
  } else {
    seurat <- Seurat::FindClusters(seurat)
  }
  return(seurat)
}



