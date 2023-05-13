#' Plot tf-idf results in reduced dimensions 
#' 
#' Plot tf-idf enrichment results in reduced dimensional space (e.g. PCA/tSNe/UMAP), 
#' Reduced dimensions can be computed based on single-cell data (e.g. RNA expression). . 
#' 
#' @inheritParams run_tfidf
#' @param size_var Point size variable in \code{obj} metadata. 
#' @param color_var Point color variable in \code{obj} metadata. 
#' @param point_alpha Point opacity.
#' @param point_palette Point palette.
#' @param density_palette Density palette.
#' @param density_adjust Density adjust (controls granularity of density plot). 
#' @param label_fill Cluster label background color.
#' @param show_plot Whether to print the plot.
#' @param background_color Plot background color.
#' @param text_color Cluster label text color.
#' @param interact Whether to make the plot interactive with \pkg{plotly}. 
#' @param ... Additional arguments to be passed to \code{ggplot2::geom_point(aes_string(...))}. 
#' 
#' @export
#' @import ggplot2
#' @examples 
#' data("pseudo_seurat")
#' res <- plot_tfidf(obj = pseudo_seurat, 
#'                   label_var = "celltype", 
#'                   cluster_var = "cluster")
plot_tfidf <- function(obj, 
                       reduction="UMAP",
                       label_var=NULL,
                       cluster_var="seurat_clusters", 
                       replace_regex="[.]|[_]|[-]",
                       terms_per_cluster=3,
                       size_var=1,
                       color_var="cluster",
                       point_alpha=.7,
                       point_palette=c(unname(pals::alphabet()),
                                       rev(unname(pals::alphabet2()) )),
                       density_palette="Purples",
                       density_adjust=.2,
                       label_fill=ggplot2::alpha(c("white"),0.7),
                       show_plot=TRUE,
                       background_color="white",
                       text_color="black",
                       interact=FALSE,
                       force_new=FALSE,
                       verbose=TRUE,
                       ...){
  # devoptera::args2vars(plot_tfidf)
  requireNamespace("ggplot2")
  
  #### Prepare clusts #### 
  r <- run_tfidf(obj=obj, 
                 reduction=reduction,
                 label_var=label_var,
                 cluster_var=cluster_var,
                 replace_regex = replace_regex,
                 terms_per_cluster=terms_per_cluster,
                 force_new=force_new,
                 return_all_results=TRUE,
                 verbose=verbose) 
   
  cluster_centers <- r$obs |>
    merge(data.frame(r$tfidf_df) |>
            # Make sure non-enriched clusters get annotated
            dplyr::mutate(word=ifelse(is.na(word),"N/A",word)),
          by="cluster") |>
    dplyr::group_by(cluster, word, .drop=FALSE) |>
    dplyr::summarise(x.mean=mean(eval(parse(text = r$dim_key[["x"]] )),
                                 na.rm=TRUE),
                     y.mean=mean(eval(parse(text = r$dim_key[["y"]] )),
                                 na.rm=TRUE),
                     size=mean(tf_idf, na.rm = TRUE)) |>
    dplyr::rename(term=word) |>
    data.table::data.table()
  
  cluster_number_centers <- cluster_centers |>
    dplyr::group_by(cluster, .drop=FALSE) |>
    dplyr::slice_head(n = 1)
  
  
  plt <- ggplot2::ggplot(data = r$obs2,
                         ggplot2::aes_string(x=r$dim_key[["x"]],
                                             y=r$dim_key[["y"]],
                                       size=size_var)) +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ggplot2::after_stat(level)),
                    adjust = density_adjust, 
                    contour = TRUE,
                    geom = "polygon",
                    contour_var = "ndensity") +
    ggplot2::scale_fill_distiller(palette=density_palette, 
                                  direction=-1) +
    ggplot2::geom_point(ggplot2::aes_string(color=color_var,
                          # shape=Type,
                          size=size_var,
                          label=label_var,),
                          # ...),
               alpha=point_alpha,
               show.legend = FALSE) +
    ggplot2::scale_color_manual(values = point_palette) +
    # geom_label(data = cluster_centers,
    #            aes(x=x.mean, y=y.mean, label=cluster, size=size, color=cluster),
    #               fill = alpha(c("white"),0.8), inherit.aes = F, show.legend = F) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = background_color),
          panel.background = ggplot2::element_rect(fill = background_color),
          legend.text = ggplot2::element_text(color=text_color),
          legend.title = ggplot2::element_text(color=text_color)
    )
  if(!is.null(label_var)){ 
    plt <- plt + 
      ggrepel::geom_label_repel(
        data = cluster_number_centers,
        ggplot2::aes(x=x.mean, y=y.mean,
                     label=cluster,
                     size=size, 
                     color=cluster),
        fill = ggplot2::alpha(c("black"),0.5),
        ## Must set this outside of aes() due to bug that causes
        ## conflict with geom_point (even when inherit.aes=F)
        size=5,
        min.segment.length = 0.1, box.padding = 3,
        inherit.aes = FALSE, max.overlaps = 30,
        show.legend = FALSE) +
      ggrepel::geom_label_repel(
        data = cluster_centers,
        ggplot2::aes(x=x.mean, y=y.mean,
                     label=term, 
                     size=size, 
                     color=cluster),
        fill = label_fill,
        ## Must set this outside of aes() due to bug that causes
        ## conflict with geom_point (even when inherit.aes=F)
        # Also, log and rescale size so that all labels are still legible
        size=scales::rescale(log10(cluster_centers$size), c(2,6)),
        min.segment.length = 0.1, box.padding = .1,
        inherit.aes = FALSE, max.overlaps = 30,
        show.legend = FALSE) 
  }
  
  if(isTRUE(interact)) {
    plt <- plotly::ggplotly(plt)
  }
  if(isTRUE(show_plot)) methods::show(plt)
  #### Return ####
  return(list(obs2=r$obs2,
              tfidf_df=r$tfidf_df,
              plot=plt))
}
