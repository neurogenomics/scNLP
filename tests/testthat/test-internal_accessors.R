test_that("get_obs_internal works", {
    testthat::skip_if_not_installed("Seurat")

    data("pseudo_seurat", package = "scNLP")

    obs <- scNLP:::get_obs_internal(pseudo_seurat, verbose = FALSE)

    # Check output is a data.frame
    testthat::expect_true(is.data.frame(obs))

    # Check it has the same number of rows as cells
    testthat::expect_equal(nrow(obs), ncol(pseudo_seurat))

    # Check rownames match cell names
    testthat::expect_equal(rownames(obs), colnames(pseudo_seurat))
})

test_that("get_obsm_internal works", {
    testthat::skip_if_not_installed("Seurat")

    data("pseudo_seurat", package = "scNLP")

    # Check UMAP exists
    testthat::skip_if_not("umap" %in% tolower(SeuratObject::Reductions(pseudo_seurat)))

    obsm <- scNLP:::get_obsm_internal(pseudo_seurat, keys = "UMAP", verbose = FALSE)

    # Check output is a list
    testthat::expect_true(is.list(obsm))

    # Check it contains the requested reduction
    testthat::expect_true(length(obsm) > 0)

    # Check the reduction has correct dimensions
    testthat::expect_equal(nrow(obsm[[1]]), ncol(pseudo_seurat))
})

test_that("set_obs_internal works", {
    testthat::skip_if_not_installed("Seurat")

    data("pseudo_seurat", package = "scNLP")

    # Get original metadata
    obs <- pseudo_seurat@meta.data

    # Add a new column
    obs$test_col <- seq_len(nrow(obs))

    # Set it back
    result <- scNLP:::set_obs_internal(pseudo_seurat, obs, verbose = FALSE)

    # Check the new column was added
    testthat::expect_true("test_col" %in% colnames(result@meta.data))
    testthat::expect_equal(result@meta.data$test_col, seq_len(ncol(result)))
})

test_that("internal accessors error on non-Seurat objects", {
    testthat::expect_error(scNLP:::get_obs_internal(data.frame()))
    testthat::expect_error(scNLP:::get_obsm_internal(data.frame()))
    testthat::expect_error(scNLP:::set_obs_internal(data.frame(), data.frame()))
})
