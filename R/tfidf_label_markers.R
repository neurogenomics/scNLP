#' TF-IDF label markers
#' 
#' Get marker genes for each cluster,
#' @inheritDotParams run_tfidf
#' 
#' @keywords internal
tfidf_label_markers <- function(obj,
                                markers,
                                ...){ 
  
  obj <- run_tfidf(obj,
                   ...)
    #### Merge markers df with cluster annotations from TF-IDF step ####
    markers_lab <- merge(
      x = markers |> dplyr::mutate(cluster=as.factor(cluster)),
      y = unique(obj@meta.data[,c("seurat_clusters","enriched_words")]), 
      all = TRUE, by.x = "cluster", by.y = "seurat_clusters")
    return(markers_lab)
}
