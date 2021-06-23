
#' Wordcloud from tf-idf results 
#' 
#' @inheritParams run_tfidf 
#' @param ... Additional parameters to pass to \code{ggplot2::ggplot(aes_string(...))}. 
#' 
#' @export
wordcloud_tfidf <- function(object,
                            label_var = "celltype", 
                            cluster_var = "cluster", 
                            terms_per_cluster=10,
                            show_plot=T,
                            ...){
  # wordcloud2: https://towardsdatascience.com/create-a-word-cloud-with-r-bde3e7422e8a
  # ggwordcloud: https://cran.r-project.org/web/packages/ggwordcloud/vignettes/ggwordcloud.html
  # object <- scNLP::pseudo_seurat
  res <- plot_tfidf(object = object, 
                    label_var = label_var, 
                    cluster_var = cluster_var, 
                    terms_per_cluster=terms_per_cluster,
                    show_plot = F) 
  dat <- res$tfidf_df
  # dat <- data.frame(res$tfidf_df, row.names = res$tfidf_df$word) %>% dplyr::select(word,freq=n)
  # plt <- wordcloud2::wordcloud2(data =  dat,
  #                               color = "random-light",
  #                               backgroundColor = "grey",
  #                               figPath = "images/brain.png")
  
  library(ggwordcloud) 
  plt <- ggplot(dat, aes_string(label = "word", size = "tf_idf", color="tf_idf", ...)) + 
    geom_text_wordcloud_area() +
    # scale_size_area(max_size = 18) +
    theme_minimal() +
    facet_wrap(~paste("cluster",cluster))
  
  if(show_plot) print(plt)
  return(list(plot=plt, 
              tfidf_df=dat))
}
