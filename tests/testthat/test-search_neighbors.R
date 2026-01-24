test_that("search_neighbors returns expected structure", {
    pseudo_seurat <- load_pseudo_seurat()

    result <- search_neighbors(
        seurat = pseudo_seurat,
        var1_search = "purkinje",
        max_neighbors = 5,
        verbose = FALSE
    )

    # Check output is a data.table/data.frame
    testthat::expect_true(is.data.frame(result))

    # Check expected columns exist
    testthat::expect_true("Var1" %in% colnames(result))
    testthat::expect_true("Var2" %in% colnames(result))
    testthat::expect_true("similarity" %in% colnames(result))
})

test_that("search_neighbors respects max_neighbors", {
    pseudo_seurat <- load_pseudo_seurat()

    result <- search_neighbors(
        seurat = pseudo_seurat,
        var1_search = "purkinje",
        max_neighbors = 3,
        verbose = FALSE
    )

    # Check that no Var1 has more than max_neighbors
    if (nrow(result) > 0) {
        neighbors_per_var1 <- table(result$Var1)
        testthat::expect_true(all(neighbors_per_var1 <= 3))
    }
})

test_that("search_neighbors filters by var2_group", {
    pseudo_seurat <- load_pseudo_seurat()

    # Skip if species column doesn't exist
    testthat::skip_if_not("species" %in% colnames(pseudo_seurat@meta.data))

    result <- search_neighbors(
        seurat = pseudo_seurat,
        var1_search = "purkinje",
        var2_group = "human",
        group_col = "species",
        max_neighbors = 5,
        verbose = FALSE
    )

    # Check that Var2_group column was added
    if (nrow(result) > 0) {
        testthat::expect_true("Var2_group" %in% colnames(result))
        # All Var2_group values should match the filter
        testthat::expect_true(all(grepl("human", result$Var2_group, ignore.case = TRUE)))
    }
})

test_that("search_neighbors adds original names when requested", {
    pseudo_seurat <- load_pseudo_seurat()

    result <- search_neighbors(
        seurat = pseudo_seurat,
        var1_search = "purkinje",
        max_neighbors = 5,
        add_original_names = TRUE,
        verbose = FALSE
    )

    # Check that original name columns were added
    if (nrow(result) > 0) {
        testthat::expect_true("Var1_id" %in% colnames(result))
        testthat::expect_true("Var2_id" %in% colnames(result))
    }
})

test_that("search_neighbors errors on no matches", {
    pseudo_seurat <- load_pseudo_seurat()

    # Search for something that definitely doesn't exist
    testthat::expect_error(
        search_neighbors(
            seurat = pseudo_seurat,
            var1_search = "xyznonexistent123",
            verbose = FALSE
        )
    )
})
