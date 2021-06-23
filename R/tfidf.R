
#' tfidf 
#' 
#' Run tf-idf on a metadata table.
#' 
#' @param clusts \code{data.frame}/\code{data.table} with the per-cell meteadata and cluster assignments. 
#' @inheritParams run_tfidf 
#' 
#' @export
tfidf <- function(clusts,
                  label_var="dataset",
                  cluster_var="seurat_clusters",
                  terms_per_cluster=1,
                  replace_regex="[.]|[_]|[-]",
                  force_new=F){
  clusts <- data.frame(clusts)
  if(!label_var %in% colnames(clusts)) stop(label_var," not found in metadata.")
  
  library(tidytext)
  data(stop_words)
  if(any(c("enriched_words","tf_idf") %in% colnames(clusts))){
      clusts <- clusts %>% dplyr::select(-c(enriched_words,tf_idf))
  }
  # From here: https://www.tidytextmining.com/tfidf.html
  clusts$cluster <- clusts[[cluster_var]]
  clusts$var <- gsub(replace_regex," ",clusts[[label_var]])
  clust_words <- clusts %>%
    unnest_tokens(word, var, drop = F) %>%
    anti_join(stop_words, keep=T) %>%
    dplyr::count(cluster, word, sort = TRUE)
  total_words <- clust_words %>%
    dplyr::group_by(cluster, .drop = F) %>%
    dplyr::summarize(total = sum(n, na.rm = T))
  total_samples <- clusts %>%
    dplyr::group_by(cluster, .drop = F) %>%
    dplyr::count(name = "samples")
  clust_words <- left_join(clust_words, total_words) %>%
    left_join(total_samples)
  clust_words <- clust_words %>%
    bind_tf_idf(word, cluster, n) %>%
    dplyr::group_by(cluster, .drop = F) %>%
    dplyr::slice_max(order_by = tf_idf, n=terms_per_cluster, with_ties = F) %>%
    subset(!word %in% c("cell","cells", stop_words))
  return(clust_words)
}


