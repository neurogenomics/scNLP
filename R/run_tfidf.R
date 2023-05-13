#' Run TF-IDF on single-cell data
#'  
#' Run Term Frequency - Inverse Document Frequency (TF-IDF) analysis on
#' samples metadata to characterise each cluster.
#' @param obj Single-cell data object.
#' @param reduction Name of the reduction to use (\emph{case insensitive}). 
#' @param label_var Which cell metadata column to input to NLP analysis. 
#' @param cluster_var Which cell metadata column to use to identify which 
#' cluster each cell is assigned to. 
#' @param replace_regex Characters by which to split \code{label_var} into terms
#'  (i.e. tokens) for NLP enrichment analysis.
#' @param terms_per_cluster The maximum number of words to return per cluster. 
#' @param force_new If NLP results are already detected the metadata, 
#' set \code{force_new=TRUE} to replace them with new results.
#' @param return_all_results Whether to return just the \code{obj} 
#' with updated metadata (\code{TRUE}), 
#' or all intermediate results (\code{FALSE}).
#' @param verbose Whether to print messages. 
#' 
#' @export
#' @examples  
#' data("pseudo_seurat")
#' obj2 <- run_tfidf(obj = pseudo_seurat, 
#'                   cluster_var = "cluster",
#'                   label_var = "celltype")  
run_tfidf <- function(obj=NULL, 
                      reduction="UMAP",
                      label_var="label",
                      cluster_var="seurat_clusters",
                      replace_regex = "[.]|[_]|[-]",
                      terms_per_cluster=3,
                      force_new=FALSE,
                      return_all_results=FALSE,
                      verbose=TRUE){
  # devoptera::args2vars(run_tfidf)
  
  #### Prepare clusts ####
  gid <- get_input_dat(obj=obj,  
                       cluster_var=cluster_var,
                       reduction=reduction,
                       verbose=verbose)
  clusts <- gid$clusts; 
  dim_key <- gid$dim_key; 
  cluster_var <- gid$cluster_var;  
  
  if(any(c("enriched_words","tf_idf") %in% colnames(clusts))){
    if(isTRUE(force_new)){
      clusts <- drop_cols(df = clusts,
                             cols = c("enriched_words","tf_idf"))
    } else {
      messager("Previous TF-IDF results detected.",
               "Use force_new=TRUE to re-run.",v=verbose)
      return(obj)
    }
  }
  tfidf_df <- tfidf(clusts = clusts,
                    label_var =  label_var,
                    cluster_var = cluster_var,
                    replace_regex = replace_regex,
                    terms_per_cluster = terms_per_cluster)
  obs2 <- merge(
    x=clusts,
    y=tfidf_df |>
      dplyr::group_by(cluster) |> 
      dplyr::summarise(enriched_words=paste(unique(word),collapse="; "),
                       tf_idf=paste(unique(tf_idf),collapse="; ")),
    all.x = TRUE,
    by.x = cluster_var, 
    by.y = "cluster",
    sort = FALSE)
  # Make sure rows are in the right order
  obs2 <- obs2[match(clusts[[label_var]], obs2[[label_var]]),]
  if(sum(obs2[[label_var]]!=clusts[[label_var]], na.rm = TRUE)>0){
    stop("Seurat sample names and tfidf samples names are not aligned!")
  }
  obs2 <- data.frame(obs2, row.names = row.names(clusts))    
  obj2 <- scKirby::set_obs(obj = obj, 
                           obs = obs2,
                           verbose = verbose)
  #### Return ####
  if(isTRUE(return_all_results)) {
    list(obj=obj2, 
         tfidf_df=tfidf_df, 
         obs=clusts,
         obs2=obs2, 
         dim_key=dim_key,
         cluster_var=cluster_var)
  } else { 
    return(obj2)
  }
}
