#' Wordcloud from tf-idf results 
#' 
#' @inheritParams run_tfidf 
#' @inheritParams plot_tfidf 
#' @param ... Additional parameters to pass to \code{ggplot2::ggplot(aes_string(...))}. 
#' 
#' @export
wordcloud_tfidf <- function(obj,
                            label_var = "celltype", 
                            cluster_var = "cluster", 
                            terms_per_cluster=10,
                            show_plot=TRUE,
                            ...){
  requireNamespace("ggwordcloud")
  requireNamespace("ggplot2")
  
  # wordcloud2: https://towardsdatascience.com/create-a-word-cloud-with-r-bde3e7422e8a
  # ggwordcloud: https://cran.r-project.org/web/packages/ggwordcloud/vignettes/ggwordcloud.html
  # obj <- scNLP::pseudo_seurat
  res <- plot_tfidf(obj = obj, 
                    label_var = label_var, 
                    cluster_var = cluster_var, 
                    terms_per_cluster=terms_per_cluster,
                    show_plot = FALSE) 
  dat <- res$tfidf_df
  # dat <- data.frame(res$tfidf_df, row.names = res$tfidf_df$word) |> dplyr::select(word,freq=n)
  # plt <- wordcloud2::wordcloud2(data =  dat,
  #                               color = "random-light",
  #                               backgroundColor = "grey",
  #                               figPath = "images/brain.png")
  

  plt <- ggplot2::ggplot(dat, ggplot2::aes_string(label = "word", 
                                size = "tf_idf",
                                color="tf_idf",
                                ...)) + 
    ggwordcloud::geom_text_wordcloud_area() +
    # scale_size_area(max_size = 18) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(~paste("cluster",cluster))
  
  if(show_plot) print(plt)
  return(list(plot=plt, 
              tfidf_df=dat))
}
