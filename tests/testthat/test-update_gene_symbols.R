test_that("update_gene_symbols correctly handles valid and invalid gene symbols", {
  test_mat <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("TP53", "BRCA1", "INVALID_GENE"),
      c("Sample1", "Sample2")
    )
  )

  result <- update_gene_symbols(test_mat, verbose = FALSE)

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)

  expect_true("TP53" %in% rownames(result))
  expect_true("BRCA1" %in% rownames(result))

  expect_false("INVALID_GENE" %in% rownames(result))
})

test_that("update_gene_symbols works with already approved symbols", {
  test_mat <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("TP53", "EGFR"),
      c("Sample1", "Sample2")
    )
  )

  result <- update_gene_symbols(test_mat, verbose = FALSE)

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)

  expect_setequal(rownames(result), c("TP53", "EGFR"))
})

test_that("update_gene_symbols throws error without rownames", {
  test_mat <- matrix(1:6, nrow = 3, ncol = 2)

  expect_error(
    update_gene_symbols(test_mat),
    "Matrix must have row names as gene identifiers"
  )
})

test_that("update_gene_symbols throws error for non-matrix input", {
  expect_error(
    update_gene_symbols(list(a = 1, b = 2)),
    "Input must be a matrix"
  )
})

test_that("update_gene_symbols handles alternate gene symbols", {
  test_mat <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("MLL", "CD33"),
      c("Sample1", "Sample2")
    )
  )

  result <- update_gene_symbols(test_mat, verbose = FALSE)

  if ("KMT2A" %in% rownames(result)) {
    expect_false("MLL" %in% rownames(result))
  }

  expect_true("CD33" %in% rownames(result))
})
