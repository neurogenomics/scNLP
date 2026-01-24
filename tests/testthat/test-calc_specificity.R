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

    # Specificity scores should sum to ~1 across columns for each gene
    row_sums <- rowSums(result)
    testthat::expect_true(all(abs(row_sums - 1) < 0.01))

    # Gene1 should be most specific to celltype_A
    testthat::expect_true(result["gene1", "celltype_A"] > result["gene1", "celltype_B"])
    testthat::expect_true(result["gene1", "celltype_A"] > result["gene1", "celltype_C"])

    # Gene3 (uniform expression) should have roughly equal specificity
    gene3_range <- max(result["gene3", ]) - min(result["gene3", ])
    testthat::expect_true(gene3_range < 0.1)
})

test_that("calc_specificity handles edge cases", {
    # Matrix with zeros
    X <- matrix(
        c(10, 0, 0,
          0, 0, 0),
        nrow = 2, byrow = TRUE
    )

    # Should not error (has small epsilon to prevent division by zero)
    result <- scNLP:::calc_specificity(X)
    testthat::expect_true(is.matrix(result))
    testthat::expect_false(any(is.nan(result)))
})
