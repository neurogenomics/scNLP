test_that("run_tfidf works with Seurat object", {
    testthat::skip_if_not_installed("Seurat")

    # Load example data
    data("pseudo_seurat", package = "scNLP")

    # Run TF-IDF
    result <- run_tfidf(
        obj = pseudo_seurat,
        cluster_var = "cluster",
        label_var = "celltype",
        verbose = FALSE
    )

    # Check output is a Seurat object
    testthat::expect_true(methods::is(result, "Seurat"))

    # Check that TF-IDF columns were added to metadata
    testthat::expect_true("enriched_words" %in% colnames(result@meta.data))
    testthat::expect_true("tf_idf" %in% colnames(result@meta.data))
})

test_that("run_tfidf respects force_new parameter", {
    testthat::skip_if_not_installed("Seurat")

    data("pseudo_seurat", package = "scNLP")

    # Run once
    result1 <- run_tfidf(
        obj = pseudo_seurat,
        cluster_var = "cluster",
        label_var = "celltype",
        verbose = FALSE
    )

    # Run again without force_new - should return early
    result2 <- run_tfidf(
        obj = result1,
        cluster_var = "cluster",
        label_var = "celltype",
        force_new = FALSE,
        verbose = FALSE
    )

    # Results should be the same (early return)
    testthat::expect_equal(
        result1@meta.data$enriched_words,
        result2@meta.data$enriched_words
    )
})
