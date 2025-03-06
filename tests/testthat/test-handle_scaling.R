data <- matrix(
  c(1, 2, 3, 4, 5, 6),
  nrow = 3,
  ncol = 2,
  dimnames = list(
    c("Gene1", "Gene2", "Gene3"),
    c("Sample1", "Sample2")
  )
)

test_that("handle_scaling with unlog = TRUE and tpm = TRUE", {
  result <- handle_scaling(
    data,
    unlog = TRUE,
    tpm = TRUE
  )
  expected <- matrix(
    c(142857.1, 285714.3, 571428.6, 142857.1, 285714.3, 571428.6),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2")
    )
  )
  expect_equal(result, expected, tolerance = 1e-1)
})

test_that("handle_scaling with unlog = FALSE and tpm = TRUE", {
  result <- handle_scaling(
    data,
    unlog = FALSE,
    tpm = TRUE
  )
  expected <- matrix(
    c(166666.7, 333333.3, 500000.0, 266666.7, 333333.3, 400000.0),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2")
    )
  )
  expect_equal(result, expected, tolerance = 1e-1)
})

test_that("handle_scaling with unlog = TRUE and tpm = FALSE", {
  result <- handle_scaling(
    data,
    unlog = TRUE,
    tpm = FALSE
  )
  expected <- matrix(
    c(2, 4, 8, 16, 32, 64),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2")
    )
  )
  expect_equal(result, expected)
})

test_that("handle_scaling with unlog = FALSE and tpm = FALSE", {
  result <- handle_scaling(
    data,
    unlog = FALSE,
    tpm = FALSE
  )
  expect_equal(result, data)
})

test_that("handle_scaling with unlog = FALSE and tpm = FALSE", {
  result <- handle_scaling(
    data,
    unlog = FALSE,
    tpm = FALSE
  )
  expect_equal(result, data)
})
