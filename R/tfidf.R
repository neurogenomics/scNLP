#' TF-IDF 
#' 
#' Run tf-idf on a metadata table. 
#' @param clusts \code{data.frame}/\code{data.table} with the 
#' per-cell metadata and cluster assignments. 
#' @inheritParams run_tfidf 
#' @inheritParams dplyr::slice_max
#' 
#' @export
tfidf <- function(clusts,
                  label_var="dataset",
                  cluster_var="seurat_clusters",
                  terms_per_cluster=1,
                  replace_regex="[.]|[_]|[-]",
                  force_new=FALSE,
                  with_ties=FALSE){
  
  requireNamespace("tidytext")
  stop_words <- tidytext::stop_words
  
  clusts <- data.frame(clusts)
  #### Check label_var ####
  if(!label_var %in% colnames(clusts)) {
    stop(label_var," not found in metadata.")
  }
  #### Drop cols ####
  clusts <- drop_cols(df = clusts,
                      cols = c("enriched_words","tf_idf"))
  # From here: https://www.tidytextmining.com/tfidf.html
  clusts$cluster <- clusts[[cluster_var]]
  clusts$var <- gsub(replace_regex," ",clusts[[label_var]])
  clust_words <- clusts |>
    tidytext::unnest_tokens(word, var, drop = FALSE) |>
    dplyr::anti_join(stop_words,by = "word") |>
    dplyr::count(cluster, word, sort = TRUE)
  total_words <- clust_words |>
    dplyr::group_by(cluster, .drop = FALSE) |>
    dplyr::summarize(total = sum(n, na.rm = TRUE))
  total_samples <- clusts |>
    dplyr::group_by(cluster, .drop = FALSE) |>
    dplyr::count(name = "samples")
  clust_words <- 
    dplyr::left_join(clust_words, total_words, by="cluster") |>
    dplyr::left_join(total_samples,by="cluster")
  clust_words <- clust_words |>
    tidytext::bind_tf_idf(word, cluster, n) |>
    dplyr::group_by(cluster, .drop = FALSE) |>
    dplyr::slice_max(order_by = tf_idf,
                     n = terms_per_cluster, 
                     with_ties = with_ties) |>
    subset(!word %in% c("cell","cells", stop_words))
  return(clust_words)
}


