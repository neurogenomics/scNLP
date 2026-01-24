test_that("seurat_pipeline works with Seurat object", {
    pseudo_seurat <- load_pseudo_seurat()

    # Get raw counts to create a fresh object
    counts <- Seurat::GetAssayData(pseudo_seurat, layer = "counts")

    # Create minimal Seurat object without reductions
    obj <- SeuratObject::CreateSeuratObject(counts = counts)

    # Run pipeline
    result <- seurat_pipeline(
        obj = obj,
        dims = 1:10,
        resolution = 0.5,
        verbose = FALSE
    )

    # Check output is a Seurat object
    testthat::expect_true(methods::is(result, "Seurat"))

    # Check that reductions were added
    testthat::expect_true("pca" %in% names(result@reductions))
    testthat::expect_true("umap" %in% names(result@reductions))

    # Check that clusters were assigned
    testthat::expect_true("seurat_clusters" %in% colnames(result@meta.data))
})

test_that("seurat_pipeline works with matrix input", {
    pseudo_seurat <- load_pseudo_seurat()

    # Get raw counts matrix
    counts <- Seurat::GetAssayData(pseudo_seurat, layer = "counts")

    # Run pipeline with matrix input
    result <- seurat_pipeline(
        obj = counts,
        dims = 1:10,
        resolution = 0.5,
        verbose = FALSE
    )

    # Check output is a Seurat object
    testthat::expect_true(methods::is(result, "Seurat"))

    # Check reductions exist
    testthat::expect_true("pca" %in% names(result@reductions))
    testthat::expect_true("umap" %in% names(result@reductions))
})

test_that("seurat_pipeline errors on invalid input", {
    testthat::expect_error(
        seurat_pipeline(obj = data.frame(a = 1:3)),
        "must be a Seurat object or a counts matrix"
    )
})
