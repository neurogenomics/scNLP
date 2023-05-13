drop_reduction_vars <- function(obs,
                                reduction="UMAP",
                                alternatives=c("umap","tsne","t-sne","t_sne","pca"),
                                verbose=TRUE){
  ## Sometimes the embeddings colnames are different from the reduction name 
  conflicting_vars <- grep(paste(c(alternatives,alternatives),collapse = "|"),
                           colnames(obs), ignore.case = TRUE, value = TRUE)
  if(length(conflicting_vars)>0){
    messager("+ Dropping",length(conflicting_vars),"conflicting obs variables:",
            paste(conflicting_vars,collapse = ", "),v=verbose)
    obs <- obs[,!colnames(obs) %in% conflicting_vars] 
  }
  return(obs)
}
