get_input_dat <- function(obj=NULL,  
                          cluster_var=NULL,
                          reduction="umap",
                          verbose=TRUE){
  # obj <- scNLP::pseudo_seurat
  
  #### Extract necessary info ####
  obs <- scKirby::get_obs(obj = obj,
                          verbose = verbose)
  obsm <- scKirby::get_obsm(obj = obj,
                            keys = reduction,
                            verbose = verbose)[[1]]  
  obs <- drop_reduction_vars(obs = obs, 
                             reduction = reduction,
                             verbose = verbose)
  clusts  <- cbind(obs, obsm)
  
  #### Infer cluster var #### 
  if(is.null(cluster_var) ||
     !cluster_var %in% colnames(clusts)){
    cluster_var <- infer_cluster_var(df = clusts)
  }
  clusts$cluster <- clusts[[cluster_var]]
  #### Create key ####
  dim_key <- setNames(colnames(obsm),c("x","y"))
  ### Drop TF-IDF vars ####
  clusts <- drop_cols(df = clusts,
                         cols = c("enriched_words","tf_idf"))
  #### Return ####
  return(list(clusts=data.frame(clusts),
              cluster_var=cluster_var,
              dim_key=dim_key))
}
