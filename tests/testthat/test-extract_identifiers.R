test_that("extract_identifiers correctly splits gene groups", {
  data <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1;Gene2;Gene3", "Gene4;Gene5", "Gene6"),
      c("Sample1", "Sample2")
    )
  )

  result <- extract_identifiers(data, separator = ";")

  expected <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene4", "Gene6"),
      c("Sample1", "Sample2")
    )
  )

  expect_equal(result, expected)
})

test_that("extract_identifiers works with custom separator", {
  data <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1|Gene2|Gene3", "Gene4|Gene5", "Gene6"),
      c("Sample1", "Sample2")
    )
  )

  result <- extract_identifiers(
    data,
    separator = "\\|"
  )

  expected <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene4", "Gene6"),
      c("Sample1", "Sample2")
    )
  )

  expect_equal(result, expected)
})

# Remove tests for NULL feature_column and invalid feature_column
# since we're using rownames in the matrix format
