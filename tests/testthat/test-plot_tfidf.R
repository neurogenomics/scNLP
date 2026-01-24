test_that("plot_tfidf returns expected structure", {
    testthat::skip_if_not_installed("ggplot2")

    pseudo_seurat <- load_pseudo_seurat()

    result <- plot_tfidf(
        obj = pseudo_seurat,
        label_var = "celltype",
        cluster_var = "cluster",
        show_plot = FALSE
    )

    # Check output is a list with expected elements
    testthat::expect_true(is.list(result))
    testthat::expect_true("data" %in% names(result))
    testthat::expect_true("tfidf_df" %in% names(result))
    testthat::expect_true("plot" %in% names(result))

    # Check plot is a ggplot object
    testthat::expect_true(methods::is(result$plot, "ggplot"))

    # Check data has expected columns
    testthat::expect_true("enriched_words" %in% colnames(result$data))
})

test_that("plot_tfidf handles custom parameters", {
    pseudo_seurat <- load_pseudo_seurat()

    result <- plot_tfidf(
        obj = pseudo_seurat,
        label_var = "celltype",
        cluster_var = "cluster",
        terms_per_cluster = 2,
        show_plot = FALSE
    )

    testthat::expect_true(is.list(result))
})
