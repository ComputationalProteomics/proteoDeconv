test_that("handle_missing_values with lowest_value imputation", {
  data <- matrix(
    c(1, NA, 3, 4, 5, NA),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2")
    )
  )
  result <- handle_missing_values(
    data,
    imputation_mode = "lowest_value"
  )
  expected <- matrix(
    c(1, 1, 3, 4, 5, 1),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2")
    )
  )
  expect_equal(result, expected)
})


test_that("handle_missing_values with knn imputation", {
  skip_if_not_installed("impute")
  data <- matrix(
    c(1, NA, 3, 4, 5, NA),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2")
    )
  )
  result <- handle_missing_values(
    data,
    imputation_mode = "knn"
  )
  expect_true(!any(is.na(result)))
})
