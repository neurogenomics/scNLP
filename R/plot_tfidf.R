
plot_tfidf <- function(object=NULL,
                       metadata=NULL,
                       embeddings=NULL,
                       reduction="UMAP",
                       label_var="label",
                       cluster_var="seurat_clusters",
                       size_var=1,
                       color_var="cluster",
                       point_alpha=.7,
                       point_palette=c(unname(pals::alphabet()),rev(unname(pals::alphabet2()) )),
                       density_palette="Purples",
                       density_adjust=.2,
                       label_fill=alpha(c("white"),0.7),
                       show_plot=T,
                       replace_regex=" ",
                       terms_per_cluster=3,
                       background_color="white",
                       text_color="black",
                       interact=F,
                       verbose=T,
                       ...){
  # object<-scNLP::pseudo_seurat;  reduction="UMAP"; label_var="label";cluster_var="seurat_clusters";size_var="genes";point_alpha=.7;density_palette="purples";density_adjust=.2; label_fill=alpha(c("white"),0.7);show_plot=T;replace_regex=" ";terms_per_cluster=3; show_plot=T; verbose=T
  
  #### Prepare input_dat #### 
  r <- run_tfidf(object=object, 
                  reduction=reduction,
                  label_var=label_var,
                  cluster_var=cluster_var,
                  replace_regex = replace_regex,
                  terms_per_cluster=terms_per_cluster,
                  force_new=F,
                  return_all_results=T,
                  verbose=verbose) 
   
  cluster_centers <- r$old_metadata %>%
    merge(data.frame(r$tfidf_df) %>%
            # Make sure non-enriched clusters get annotated
            dplyr::mutate(word=ifelse(is.na(word),"N/A",word)),
          by="cluster") %>%
    dplyr::group_by(cluster, word, .drop=F) %>%
    dplyr::summarise(x.mean=mean(eval(parse(text = r$dim_key[["x"]] )), na.rm=T),
                     y.mean=mean(eval(parse(text = r$dim_key[["y"]] )), na.rm=T),
                     size=mean(tf_idf, na.rm = T)) %>%
    dplyr::rename(term=word) %>%
    data.table::data.table()
  
  cluster_number_centers <- cluster_centers %>%
    group_by(cluster, .drop=F) %>%
    dplyr::slice_head(n = 1)
  
  library(ggplot2)
  umap_tfidf_plot <- ggplot(data = r$new_metadata,
                            aes_string(x=r$dim_key[["x"]], y=r$dim_key[["y"]], size=size_var)) +
    stat_density_2d(aes(fill = ..level..),
                    adjust = density_adjust, contour = T,
                    geom = "polygon",contour_var = "ndensity") +
    scale_fill_distiller(palette=density_palette, direction=-1) +
    geom_point(aes_string(color=color_var,
                          # shape=Type,
                          size=size_var,
                          label=label_var,
                          ...),
               alpha=point_alpha, show.legend = F) +
    scale_color_manual(values = point_palette) +
    # geom_label(data = cluster_centers,
    #            aes(x=x.mean, y=y.mean, label=cluster, size=size, color=cluster),
    #               fill = alpha(c("white"),0.8), inherit.aes = F, show.legend = F) +
    ggrepel::geom_label_repel(data = cluster_number_centers,
                              aes(x=x.mean, y=y.mean, label=cluster, size=size, color=cluster),
                              fill = alpha(c("black"),0.5),
                              ## Must set this outside of aes() due to bug that causes
                              ## conflict with geom_point (even when inherit.aes=F)
                              size=5,
                              min.segment.length = 0.1, box.padding = 3,
                              inherit.aes = F, max.overlaps = 30, show.legend = F) +
    ggrepel::geom_label_repel(data = cluster_centers,
                              aes(x=x.mean, y=y.mean, label=term, size=size, color=cluster),
                              fill = label_fill,
                              ## Must set this outside of aes() due to bug that causes
                              ## conflict with geom_point (even when inherit.aes=F)
                              # Also, log and rescale size so that all labels are still legible
                              size=scales::rescale(log10(cluster_centers$size), c(2,6)),
                              min.segment.length = 0.1, box.padding = .1,
                              inherit.aes = F, max.overlaps = 30, show.legend = F) +
    theme_void() +
    theme(plot.background = element_rect(fill = background_color),
          panel.background = element_rect(fill = background_color),
          legend.text = element_text(color=text_color),
          legend.title = element_text(color=text_color)
    )
  
  if(interact) umap_tfidf_plot <- plotly::ggplotly(umap_tfidf_plot)
  if(show_plot) print(umap_tfidf_plot)
  
  return(list(new_metadata=r$new_metadata,
              tfidf_df=r$tfidf_df,
              plot=umap_tfidf_plot))
}