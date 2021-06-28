
drop_reduction_vars <- function(metadata,
                                reduction="UMAP",
                                alternatives=c("umap","tsne","t-sne","t_sne","pca"),
                                verbose=T){
  ## Sometimes the embeddings colnames are different from the reduction name 
  conflicting_vars <- grep(paste(c(alternatives,alternatives),collapse = "|"), colnames(metadata), ignore.case = T, value = T)
  if(length(conflicting_vars)>0){
    printer("+ Dropping",length(conflicting_vars),"conflicting metadata variables:",paste(conflicting_vars,collapse = ", "),v=verbose)
    metadata <- metadata[,!colnames(metadata) %in% conflicting_vars] 
  }
  return(metadata)
}
