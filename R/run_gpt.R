#' Run GPT on single-cell data
#' 
#' Ask chatGPT to summarise each cluster based on the samples metadata.
#' @inheritParams run_tfidf
#' 
#' @export
#' @examples  
#' data("pseudo_seurat")
#' obj2 <- run_gpt(obj = pseudo_seurat, 
#'                 cluster_var = "cluster",
#'                 label_var = "celltype")  
run_gpt <- function(obj=NULL, 
                    reduction="UMAP",
                    label_var="label",
                    cluster_var="seurat_clusters", 
                    terms_per_cluster=3,
                    force_new=FALSE,
                    return_all_results=FALSE,
                    verbose=TRUE){
  # devoptera::args2vars(run_gpt)
  
  #### Prepare clusts ####
  gid <- get_input_dat(obj=obj,  
                       cluster_var=cluster_var,
                       reduction=reduction,
                       verbose=verbose)
  clusts <- gid$clusts; 
  dim_key <- gid$dim_key; 
  cluster_var <- gid$cluster_var;  
  
  if(any(c("gpt_summary") %in% colnames(clusts))){
    if(isTRUE(force_new)){
      clusts <- drop_cols(df = clusts,
                          cols =c("gpt_summary","gpt"))
    } else {
      messager("Previous GPT results detected.",
               "Use force_new=TRUE to re-run.",v=verbose)
      return(obj)
    }
  } 
  gpt_df <- gpt(clusts = clusts,
                label_var =  label_var,
                cluster_var = cluster_var, 
                terms_per_cluster = terms_per_cluster, 
                verbose = verbose)
  obs2 <- merge(
    x = clusts,
    y = gpt_df[,c("cluster","gpt_summary")],
    all.x = TRUE,
    by.x = cluster_var, 
    by.y = "cluster",
    sort = FALSE)
  # Make sure rows are in the right order
  obs2 <- obs2[match(clusts[[label_var]], obs2[[label_var]]),]
  if(sum(obs2[[label_var]]!=clusts[[label_var]], na.rm = TRUE)>0){
    stop("Seurat sample names and GPT sample names are not aligned!")
  }
  obs2 <- data.frame(obs2, row.names = row.names(clusts))    
  obj2 <- scKirby::set_obs(obj = obj, 
                           obs = obs2,
                           verbose = verbose)
  #### Return ####
  if(isTRUE(return_all_results)) {
    list(obj=obj2, 
         gpt_df=gpt_df, 
         obs=clusts,
         obs2=obs2, 
         dim_key=dim_key,
         cluster_var=cluster_var)
  } else { 
    return(obj2)
  }
}
