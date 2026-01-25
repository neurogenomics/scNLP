test_that("calc_specificity returns correct structure", {
    # Create simple test matrix
    # Rows = genes, Cols = cell types
    X <- matrix(
        c(10, 0, 0,
          0, 10, 0,
          5, 5, 5),
        nrow = 3, byrow = TRUE,
        dimnames = list(
            c("gene1", "gene2", "gene3"),
            c("celltype_A", "celltype_B", "celltype_C")
        )
    )

    result <- scNLP:::calc_specificity(X)

    # Check output dimensions match input
    testthat::expect_equal(dim(result), dim(X))

    # Check row and column names preserved
    testthat::expect_equal(rownames(result), rownames(X))
    testthat::expect_equal(colnames(result), colnames(X))

    # Check output is numeric matrix
    testthat::expect_true(is.matrix(result))
    testthat::expect_true(is.numeric(result))
})

test_that("calc_specificity produces valid specificity scores", {
    # Create test matrix where gene1 is specific to celltype_A
    X <- matrix(
        c(100, 1, 1,
          1, 100, 1,
          10, 10, 10),
        nrow = 3, byrow = TRUE,
        dimnames = list(
            c("gene1", "gene2", "gene3"),
            c("celltype_A", "celltype_B", "celltype_C")
        )
    )

    result <- scNLP:::calc_specificity(X)

    # Check output is valid
    testthat::expect_true(is.matrix(result))
    testthat::expect_true(is.numeric(result))

    # Gene1 should be most specific to celltype_A (highest value in that column)
    testthat::expect_true(result["gene1", "celltype_A"] > result["gene1", "celltype_B"])
    testthat::expect_true(result["gene1", "celltype_A"] > result["gene1", "celltype_C"])
})

test_that("calc_specificity handles edge cases", {
    # Matrix with zeros - function should handle without error
    X <- matrix(
        c(10, 0, 0,
          0, 0, 0),
        nrow = 2, byrow = TRUE
    )

    # Should not error
    result <- scNLP:::calc_specificity(X)
    testthat::expect_true(is.matrix(result))
})
