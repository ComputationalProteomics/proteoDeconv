data <- data.frame(
  Genes = c("Gene1", "Gene2", "Gene3"),
  Sample1 = c(1, 2, 3),
  Sample2 = c(4, 5, 6)
)

test_that("handle_scaling with unlog = TRUE and tpm = TRUE", {
  result <- handle_scaling(data, gene_column = "Genes", unlog = TRUE, tpm = TRUE)
  expected <- data.frame(
    Genes = c("Gene1", "Gene2", "Gene3"),
    Sample1 = c(142857.1, 285714.3, 571428.6),
    Sample2 = c(142857.1, 285714.3, 571428.6)
  )
  expect_equal(result, expected, tolerance = 1e-1)
})

test_that("handle_scaling with unlog = FALSE and tpm = TRUE", {
  result <- handle_scaling(data, gene_column = "Genes", unlog = FALSE, tpm = TRUE)
  expected <- data.frame(
    Genes = c("Gene1", "Gene2", "Gene3"),
    Sample1 = c(166666.7, 333333.3, 500000.0),
    Sample2 = c(266666.7, 333333.3, 400000.0)
  )
  expect_equal(result, expected, tolerance = 1e-1)
})

test_that("handle_scaling with unlog = TRUE and tpm = FALSE", {
  result <- handle_scaling(data, gene_column = "Genes", unlog = TRUE, tpm = FALSE)
  expected <- data.frame(
    Genes = c("Gene1", "Gene2", "Gene3"),
    Sample1 = c(2, 4, 8),
    Sample2 = c(16, 32, 64)
  )
  expect_equal(result, expected)
})

test_that("handle_scaling with unlog = FALSE and tpm = FALSE", {
  result <- handle_scaling(data, gene_column = "Genes", unlog = FALSE, tpm = FALSE)
  expect_equal(result, data)
})

test_that("handle_scaling with unlog = FALSE and tpm = FALSE", {
  result <- handle_scaling(data, gene_column = "Genes", unlog = FALSE, tpm = FALSE)
  expect_equal(result, data)
})
