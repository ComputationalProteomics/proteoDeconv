test_that("simulate_data validates input correctly", {
  data <- matrix(
    c(150, 200, 250, 300, 350, 400),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2")
    )
  )
  cell_types <- c("TypeA")
  
  expect_error(
    simulate_data(data, cell_types),
    "Length of cell_types must match the number of samples"
  )
})

test_that("simulate_data produces expected output structure", {
  data <- matrix(
    c(100, 200, 300, 400, 500, 600, 700, 800, 900),
    nrow = 3,
    ncol = 3,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2", "Sample3")
    )
  )
  cell_types <- c("TypeA", "TypeB", "TypeC")
  
  result <- simulate_data(data, cell_types, seed = 42, nsamples = 5)
  
  expect_true(is.list(result))
  expect_true(all(c("simulated_data", "cell_fractions") %in% names(result)))
  
  expect_true(is.matrix(result$simulated_data))
  expect_equal(nrow(result$simulated_data), 3)
  expect_equal(ncol(result$simulated_data), 5)
  expect_equal(rownames(result$simulated_data), c("Gene1", "Gene2", "Gene3"))
  
  expect_true(all(cell_types %in% colnames(result$cell_fractions)))
  expect_equal(nrow(result$cell_fractions), 5)
})

test_that("simulate_data creates requested number of samples", {
  data <- matrix(
    c(1000, 2000, 3000, 4000),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("Sample1", "Sample2")
    )
  )
  cell_types <- c("TypeA", "TypeB")
  
  result <- simulate_data(data, cell_types, nsamples = 8, seed = 123)
  
  expect_equal(nrow(result$cell_fractions), 8)
  expect_equal(ncol(result$simulated_data), 8)
})

test_that("simulate_data results are reproducible with seed", {
  data <- matrix(
    c(1000, 2000, 3000, 4000),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("Sample1", "Sample2")
    )
  )
  cell_types <- c("TypeA", "TypeB")
  
  result1 <- simulate_data(data, cell_types, seed = 42, nsamples = 3)
  result2 <- simulate_data(data, cell_types, seed = 42, nsamples = 3)
  
  expect_equal(result1$simulated_data, result2$simulated_data)
  expect_equal(result1$cell_fractions, result2$cell_fractions)
})

test_that("simulate_data respects whitelist and blacklist", {
  data <- matrix(
    c(100, 200, 300, 400, 500, 600, 700, 800),
    nrow = 2,
    ncol = 4,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("Sample1", "Sample2", "Sample3", "Sample4")
    )
  )
  cell_types <- c("TypeA", "TypeB", "TypeC", "TypeD")
  
  result_white <- simulate_data(
    data,
    cell_types,
    whitelist = c("TypeA", "TypeB"),
    seed = 42,
    nsamples = 3
  )
  
  expect_true(all(
    c("TypeA", "TypeB") %in% colnames(result_white$cell_fractions)
  ))
  expect_false(any(
    c("TypeC", "TypeD") %in% colnames(result_white$cell_fractions)
  ))
})