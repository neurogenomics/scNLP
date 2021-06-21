

infer_cluster_var <- function(df,
                              verbose=T){
  # df <- seurat@meta.data 
  printer("+ Inferring cluster_var..",v=verbose)
  ### Attempt 1
  cluster_var <- grep("seurat_clusters",colnames(df), ignore.case = T, value = T)[1]
  if(length(cluster_var)>0) return(cluster_var)
  ### Attempt 2
  cluster_var <- grep("cluster",colnames(df), ignore.case = T, value = T)[1]
  if(length(cluster_var)>0) return(cluster_var)
  ### Attempt 3
  cluster_var <- grep("clust",colnames(df), ignore.case = T, value = T)[1]
  if(length(cluster_var)>0) return(cluster_var) 
}
