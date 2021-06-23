
#' Run standardized \pkg{Seurat} pipeline
#' 
#' Run \pkg{Seurat} pipeline on \pkg{Seurat} object 
#' or raw \code{counts} and \code{meta.data}. 
#' 
#' Automatically performs 
#' \describe{
#' \item{\code{FindVariableFeatures}}{Variable feature selection}
#' \item{\code{NormalizeData}}{Data normalization}
#' \item{\code{RunPCA}}{PCA}
#' \item{\code{RunUMAP}}{UMAP}
#' \item{\code{FindNeighbors}}{K-nearest neighbors}
#' \item{\code{FindClusters}}{Clustering}
#' }
#' 
#' @export
seurat_pipeline <- function(seurat_obj=NULL,
                            counts=NULL, meta.data=NULL,
                            nfeatures=2000,
                            vars.to.regress=NULL,
                            dims=1:50,
                            add_specificity=F,
                            assay_name="RNA",
                            default_assay=NULL,
                            n.components=2L,
                            log_norm=F,
                            parallelize=T,
                            seed=2020){
  set.seed(seed)
  library(Seurat)
  library(future)
  if(parallelize) plan(strategy = "multicore", workers = future::availableCores()-2)
  if(is.null(seurat_obj)){
    seurat <- Seurat::CreateSeuratObject(counts = counts,  
                                         meta.data = meta.data, 
                                         assay = assay_name 
    )
  } else { seurat <- seurat_obj}
  
  if(add_specificity){
    print("+ Calculating specificity as a new assay...") 
    spec1 <- calc_specificity(as.matrix(Seurat::GetAssayData(seurat)))
    seurat[['specificity']] <- Seurat::CreateAssayObject(counts = spec1) 
  }
  
  ### Set assay ###
  if(!is.null(default_assay)) Seurat::DefaultAssay(seurat) <- default_assay 
  #### Select variable features ####
  seurat <- Seurat::FindVariableFeatures(seurat, nfeatures=nfeatures)
  seurat <- Seurat::NormalizeData(seurat)
  if(log_norm) logged <- Seurat::LogNormalize(seurat)
  seurat <- Seurat::ScaleData(seurat,  vars.to.regress = vars.to.regress)
  
  
  seurat <- Seurat::RunPCA(seurat)
  seurat <- Seurat::FindNeighbors(seurat)
  seurat <- Seurat::RunUMAP(seurat, dims=dims, n.components = n.components)
  ### Recompute neighbors graph based on UMAP to get clusters that best reflect UMAP space.
  #### For this same reason, only cluster in two dimensions, bc this is the view we most often use
  seurat <- Seurat::FindNeighbors(seurat, reduction="umap", dims=1:2)
  seurat <- Seurat::FindClusters(seurat, reduction="umap")
  return(seurat)
}



