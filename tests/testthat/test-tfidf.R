test_that("tfidf works", {
    # Create simple test data
    clusts <- data.frame(
        cluster = c(1, 1, 1, 2, 2, 2),
        celltype = c("Neuron_A", "Neuron_B", "Neuron_C",
                     "Astrocyte_A", "Astrocyte_B", "Glia_C"),
        row.names = paste0("cell", 1:6)
    )

    result <- tfidf(
        clusts = clusts,
        label_var = "celltype",
        cluster_var = "cluster",
        terms_per_cluster = 2
    )

    # Check output structure
    testthat::expect_true(is.data.frame(result))
    testthat::expect_true("word" %in% colnames(result))
    testthat::expect_true("tf_idf" %in% colnames(result))
    testthat::expect_true("cluster" %in% colnames(result))

    # Check that we get results for both clusters
    testthat::expect_true(all(c(1, 2) %in% result$cluster))
})

test_that("tfidf handles edge cases", {
    # Single cluster
    clusts <- data.frame(
        cluster = c(1, 1, 1),
        celltype = c("TypeA", "TypeB", "TypeC"),
        row.names = paste0("cell", 1:3)
    )

    # Should still work with single cluster
    result <- tfidf(
        clusts = clusts,
        label_var = "celltype",
        cluster_var = "cluster"
    )

    testthat::expect_true(is.data.frame(result))
})
