test_that("wordcloud_tfidf returns expected structure", {
    testthat::skip_if_not_installed("ggwordcloud")
    testthat::skip_if_not_installed("ggplot2")

    pseudo_seurat <- load_pseudo_seurat()

    result <- wordcloud_tfidf(
        obj = pseudo_seurat,
        label_var = "celltype",
        cluster_var = "cluster",
        terms_per_cluster = 5,
        show_plot = FALSE
    )

    # Check output is a list with expected elements
    testthat::expect_true(is.list(result))
    testthat::expect_true("plot" %in% names(result))
    testthat::expect_true("tfidf_df" %in% names(result))

    # Check plot is a ggplot object
    testthat::expect_true(methods::is(result$plot, "ggplot"))

    # Check tfidf_df has expected columns
    testthat::expect_true("word" %in% colnames(result$tfidf_df))
    testthat::expect_true("tf_idf" %in% colnames(result$tfidf_df))
    testthat::expect_true("cluster" %in% colnames(result$tfidf_df))
})

test_that("wordcloud_tfidf respects terms_per_cluster", {
    testthat::skip_if_not_installed("ggwordcloud")

    pseudo_seurat <- load_pseudo_seurat()

    result <- wordcloud_tfidf(
        obj = pseudo_seurat,
        label_var = "celltype",
        cluster_var = "cluster",
        terms_per_cluster = 3,
        show_plot = FALSE
    )

    # Check that we don't have more than 3 terms per cluster
    terms_per_clust <- table(result$tfidf_df$cluster)
    testthat::expect_true(all(terms_per_clust <= 3))
})
