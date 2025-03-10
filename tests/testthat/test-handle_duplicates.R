test_that("duplication mode=slice works", {
  data <- matrix(
    c(1, 2, 1, 2, 3, 2, 3, 3, 3),
    nrow = 3,
    ncol = 3,
    dimnames = list(
      c("A", "A", "C"),
      c("sample1", "sample2", "sample3")
    )
  )
  data_expected <- matrix(
    c(2, 1, 3, 2, 3, 3),
    nrow = 2,
    ncol = 3,
    dimnames = list(
      c("A", "C"),
      c("sample1", "sample2", "sample3")
    )
  )
  expect_equal(
    handle_duplicates(data, duplicate_mode = "slice"),
    data_expected
  )
})
