#' GPT
#' 
#' Run GPT on a metadata table. 
#' @param clusts \code{data.frame}/\code{data.table} with the 
#' per-cell metadata and cluster assignments. 
#' @inheritParams run_tfidf 
#' @inheritParams dplyr::slice_max
#' 
#' @export
#' @import data.table
gpt <- function(clusts,
                label_var="dataset",
                cluster_var="seurat_clusters",
                terms_per_cluster=5, 
                force_new=FALSE, 
                verbose=TRUE){
  
  requireNamespace("gptstudio")  
  clusts <- data.frame(clusts)
  #### Check label_var ####
  if(!label_var %in% colnames(clusts)) {
    stop(label_var," not found in metadata.")
  }
  #### Drop cols ####
  clusts <- drop_cols(df = clusts,
                      cols = c("gpt_summary"))
  clusts$cluster <- clusts[[cluster_var]]
  gpt_df <- lapply(unique(clusts$cluster), function(x){
    messager("GPT: summarising cluster",x,v=verbose)
    prompt <- paste(
      "Summarise the following text in less than",terms_per_cluster,"words:",
      paste(unique(subset(clusts,cluster==x)[[label_var]]),collapse = ". ")
    )  
    out <- gptstudio::openai_create_chat_completion(prompt = prompt)
    data.table::data.table( 
      prompt = prompt,
      cluster = x,
      data.table::as.data.table(out$usage),
      out$choices
    )
  }) |> data.table::rbindlist()   
  data.table::setnames(gpt_df,"message.content","gpt_summary")
  return(gpt_df)
}


