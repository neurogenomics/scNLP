
#' Get marker genes for each cluster 
tfidf_label_markers <- function(seurat,
                                markers,
                                label_var=NULL){
    if(!all(c("seurat_clusters","enriched_words") %in% colnames(seurat@meta.data))){ 
        if(is.null(label_var)) stop("Provide `label_var`")
        seurat <- seurat_tfidf(seurat,
                               label_var=label_var,
                               cluster_var="seurat_clusters", 
                               terms_per_cluster=3,
                               force_new = T)
    }
    #### Merge markers df with cluster annotations from TF-IDF step ####
    markers_lab <- merge(markers |> dplyr::mutate(cluster=as.factor(cluster)),
                         unique(seurat@meta.data[,c("seurat_clusters","enriched_words")]), 
                         all = T, by.x = "cluster", by.y = "seurat_clusters")
    return(markers_lab)
}
