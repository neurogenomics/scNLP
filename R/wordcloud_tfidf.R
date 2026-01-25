#' Wordcloud from tf-idf results
#'
#' @inheritParams run_tfidf
#' @inheritParams plot_tfidf
#' @param ... Additional parameters passed to \code{ggplot2::ggplot(aes(...))}.
#'
#' @returns A list containing:
#' \describe{
#'   \item{plot}{The wordcloud ggplot object.}
#'   \item{tfidf_df}{The TF-IDF results data.frame.}
#' }
#' @export
#' @examples
#' data("pseudo_seurat")
#' if (requireNamespace("ggwordcloud", quietly = TRUE)) {
#'     wordcloud_res <- wordcloud_tfidf(obj = pseudo_seurat,
#'                                      label_var = "celltype",
#'                                      cluster_var = "cluster")
#' }
wordcloud_tfidf <- function(obj,
                            label_var = "celltype",
                            cluster_var = "cluster",
                            terms_per_cluster = 10,
                            show_plot = TRUE,
                            ...) {
    requireNamespace("ggwordcloud")
    requireNamespace("ggplot2")

    res <- plot_tfidf(obj = obj,
                      label_var = label_var,
                      cluster_var = cluster_var,
                      terms_per_cluster = terms_per_cluster,
                      show_plot = FALSE)
    dat <- res$tfidf_df

    plt <- ggplot2::ggplot(
        dat,
        ggplot2::aes(
            label = .data$word,
            size = .data$tf_idf,
            color = .data$tf_idf,
            ...
        )
    ) +
        ggwordcloud::geom_text_wordcloud_area() +
        ggplot2::theme_minimal() +
        ggplot2::facet_wrap(~ paste("cluster", cluster))

    if (show_plot) methods::show(plt)
    return(list(plot = plt,
                tfidf_df = dat))
}
