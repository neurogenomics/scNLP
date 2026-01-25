test_that("FindVariableFeatures_split works", {
    pseudo_seurat <- load_pseudo_seurat()

    # Check that 'batch' column exists
    testthat::skip_if_not("batch" %in% colnames(pseudo_seurat@meta.data),
                         message = "batch column not in pseudo_seurat")

    var_features <- FindVariableFeatures_split(
        seurat = pseudo_seurat,
        split.by = "batch",
        nfeatures = 100,
        nfeatures_max = 50
    )

    # Check output is a character vector
    testthat::expect_true(is.character(var_features))

    # Check we get the expected number of features (or less)
    testthat::expect_lte(length(var_features), 50)

    # Check features are valid gene names from the object
    testthat::expect_true(all(var_features %in% rownames(pseudo_seurat)))
})

test_that("FindVariableFeatures_split returns nested list when requested", {
    pseudo_seurat <- load_pseudo_seurat()
    testthat::skip_if_not("batch" %in% colnames(pseudo_seurat@meta.data))

    var_features <- FindVariableFeatures_split(
        seurat = pseudo_seurat,
        split.by = "batch",
        nfeatures = 100,
        return_nested = TRUE
    )

    # Check output is a list
    testthat::expect_true(is.list(var_features))

    # Check names match batches
    batches <- unique(pseudo_seurat@meta.data$batch)
    testthat::expect_true(all(names(var_features) %in% batches))
})
