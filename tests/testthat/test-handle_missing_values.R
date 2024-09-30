test_that("handle_missing_values with lowest_value imputation", {
  data <- data.frame(
    Genes = c("Gene1", "Gene2", "Gene3"),
    Sample1 = c(1, NA, 3),
    Sample2 = c(4, 5, NA)
  )
  result <- handle_missing_values(data, gene_column = "Genes", imputation_mode = "lowest_value")
  expected <- data.frame(
    Genes = c("Gene1", "Gene2", "Gene3"),
    Sample1 = c(1, 1, 3),
    Sample2 = c(4, 5, 1)
  )
  expect_equal(result, expected)
})


test_that("handle_missing_values with knn imputation", {
  data <- data.frame(
    Genes = c("Gene1", "Gene2", "Gene3"),
    Sample1 = c(1, NA, 3),
    Sample2 = c(4, 5, NA)
  )
  result <- handle_missing_values(data, gene_column = "Genes", imputation_mode = "knn")
  expect_true(!any(is.na(result)))
})

