infer_cluster_var <- function(df,
                              verbose=TRUE){
  # df <- seurat@meta.data 
  
  messager("+ Inferring cluster_var..",v=verbose)
  ### Attempt 1
  cluster_var <- grep("seurat_clusters",colnames(df),
                      ignore.case = TRUE, value = TRUE)[1]
  if(length(cluster_var)>0) return(cluster_var)
  ### Attempt 2
  cluster_var <- grep("cluster",colnames(df),
                      ignore.case = TRUE, value = TRUE)[1]
  if(length(cluster_var)>0) return(cluster_var)
  ### Attempt 3
  cluster_var <- grep("clust",colnames(df),
                      ignore.case = TRUE, value = TRUE)[1]
  if(length(cluster_var)>0) return(cluster_var) 
}
