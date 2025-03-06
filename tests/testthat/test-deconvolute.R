test_that("deconvolute loads and accesses test data correctly", {
  pure_data_file <- system.file(
    "extdata",
    "pure_samples_matrix.rds",
    package = "proteoDeconv"
  )
  expect_false(pure_data_file == "", "cd8_mono.rds not found")

  mix_data_file <- system.file(
    "extdata",
    "mixed_samples_matrix.rds",
    package = "proteoDeconv"
  )
  expect_false(mix_data_file == "", "mixed_samples_matrix.rds not found")

  signature_file <- system.file(
    "extdata",
    "cd8t_mono_signature_matrix.rds",
    package = "proteoDeconv"
  )
  expect_false(signature_file == "", "cd8t_mono_signature_matrix.rds not found")

  pure_samples <- readRDS(pure_data_file)
  if (is.data.frame(pure_samples)) {
    gene_col <- which(sapply(pure_samples, is.character))[1]
    genes <- pure_samples[[gene_col]]
    mat_data <- as.matrix(pure_samples[, -gene_col])
    rownames(mat_data) <- genes
    pure_samples <- mat_data
  }
  expect_true(is.matrix(pure_samples), "pure_samples should be a matrix")

  mixed_samples <- readRDS(mix_data_file)
  if (is.data.frame(mixed_samples)) {
    gene_col <- which(sapply(mixed_samples, is.character))[1]
    genes <- mixed_samples[[gene_col]]
    mat_data <- as.matrix(mixed_samples[, -gene_col])
    rownames(mat_data) <- genes
    mixed_samples <- mat_data
  }
  expect_true(is.matrix(mixed_samples), "mixed_samples should be a matrix")

  signature <- readRDS(signature_file)
  if (is.data.frame(signature)) {
    gene_col <- which(sapply(signature, is.character))[1]
    genes <- signature[[gene_col]]
    mat_data <- as.matrix(signature[, -gene_col])
    rownames(mat_data) <- genes
    signature <- mat_data
  }
  expect_true(is.matrix(signature), "signature should be a matrix")
})

test_that("deconvolute preprocessing functions work correctly", {
  mix_data_file <- system.file(
    "extdata",
    "mixed_samples_matrix.rds",
    package = "proteoDeconv"
  )
  mixed_samples <- readRDS(mix_data_file)

  if (is.data.frame(mixed_samples)) {
    gene_col <- which(sapply(mixed_samples, is.character))[1]
    genes <- mixed_samples[[gene_col]]
    mat_data <- as.matrix(mixed_samples[, -gene_col])
    rownames(mat_data) <- genes
    mixed_samples <- mat_data
  }

  expect_error(
    preprocessed_samples <- mixed_samples |>
      extract_identifiers() |>
      update_gene_symbols(verbose = FALSE) |>
      handle_missing_values(imputation_mode = "lowest_value") |>
      handle_duplicates(duplicate_mode = "slice") |>
      handle_scaling(unlog = FALSE, tpm = TRUE),
    regexp = NA
  )

  if (exists("preprocessed_samples")) {
    expect_true(is.matrix(preprocessed_samples))
  }
})

test_that("deconvolute accepts signature as matrix", {
  bulk_matrix <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Sample1", "Sample2")
    )
  )

  signature <- matrix(
    c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("CellType1", "CellType2")
    )
  )

  expect_error(
    deconvolute(
      algorithm = "mock",
      data = bulk_matrix,
      signature = signature
    ),
    "Unknown algorithm"
  )
})

test_that("deconvolute EPIC works correctly", {
  skip_if_not_installed("EPIC")

  mix_data_file <- system.file(
    "extdata",
    "mixed_samples_matrix.rds",
    package = "proteoDeconv"
  )
  signature_file <- system.file(
    "extdata",
    "cd8t_mono_signature_matrix.rds",
    package = "proteoDeconv"
  )

  mixed_samples <- readRDS(mix_data_file)
  signature_matrix <- readRDS(signature_file)

  if (is.data.frame(mixed_samples)) {
    gene_col <- which(sapply(mixed_samples, is.character))[1]
    genes <- mixed_samples[[gene_col]]
    mat_data <- as.matrix(mixed_samples[, -gene_col])
    rownames(mat_data) <- genes
    mixed_samples <- mat_data
  }

  if (is.data.frame(signature_matrix)) {
    gene_col <- which(sapply(signature_matrix, is.character))[1]
    genes <- signature_matrix[[gene_col]]
    mat_data <- as.matrix(signature_matrix[, -gene_col])
    rownames(mat_data) <- genes
    signature_matrix <- mat_data
  }

  preprocessed_samples <- mixed_samples |>
    extract_identifiers() |>
    update_gene_symbols(verbose = FALSE) |>
    handle_missing_values(imputation_mode = "lowest_value") |>
    handle_duplicates(duplicate_mode = "slice") |>
    handle_scaling(unlog = FALSE, tpm = TRUE)

  results <- deconvolute(
    algorithm = "epic",
    data = preprocessed_samples,
    signature = signature_matrix
  )

  expect_true(is.matrix(results))
  expect_true(ncol(results) >= 1)
  expect_true(nrow(results) >= 1)

  expected_cell_types <- colnames(signature_matrix)
  cell_types_in_results <- colnames(results)
  expect_true(any(expected_cell_types %in% cell_types_in_results))
})

test_that("deconvolute BayesDebulk works correctly", {
  skip_if_not_installed("BayesDeBulk")

  mix_data_file <- system.file(
    "extdata",
    "mixed_samples_matrix.rds",
    package = "proteoDeconv"
  )

  mixed_samples <- readRDS(mix_data_file)
  if (is.data.frame(mixed_samples)) {
    gene_col <- which(sapply(mixed_samples, is.character))[1]
    genes <- mixed_samples[[gene_col]]
    mat_data <- as.matrix(mixed_samples[, -gene_col])
    rownames(mat_data) <- genes
    mixed_samples <- mat_data
  }

  signature_matrix <- readRDS(system.file(
    "extdata",
    "cd8t_mono_signature_matrix.rds",
    package = "proteoDeconv"
  ))

  if (is.data.frame(signature_matrix)) {
    gene_col <- which(sapply(signature_matrix, is.character))[1]
    genes <- signature_matrix[[gene_col]]
    mat_data <- as.matrix(signature_matrix[, -gene_col])
    rownames(mat_data) <- genes
    signature_matrix <- mat_data
  }
  signature_matrix <- handle_duplicates(
    signature_matrix,
    duplicate_mode = "slice"
  )

  preprocessed_samples <- mixed_samples |>
    extract_identifiers() |>
    update_gene_symbols(verbose = FALSE) |>
    handle_missing_values(imputation_mode = "lowest_value") |>
    handle_duplicates(duplicate_mode = "slice") |>
    handle_scaling(unlog = FALSE, tpm = TRUE)

  results <- deconvolute(
    algorithm = "bayesdebulk",
    data = preprocessed_samples,
    signature = signature_matrix
  )

  expect_true(is.matrix(results))
  expect_true(ncol(results) >= 1)
  expect_true(nrow(results) >= 1)

  expect_false(is.null(colnames(results)))
})

test_that("deconvolute CIBERSORTx fails appropriately when environment variables are not set", {
  bulk_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("Sample1", "Sample2")
    )
  )

  signature <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("CellType1", "CellType2")
    )
  )

  old_email <- Sys.getenv("CIBERSORTX_EMAIL")
  old_token <- Sys.getenv("CIBERSORTX_TOKEN")

  Sys.unsetenv("CIBERSORTX_EMAIL")
  Sys.unsetenv("CIBERSORTX_TOKEN")

  expect_error(
    deconvolute(
      algorithm = "cibersortx",
      data = bulk_matrix,
      signature = signature
    ),
    "CIBERSORTX_EMAIL and CIBERSORTX_TOKEN environment variables must be set"
  )

  if (old_email != "") Sys.setenv(CIBERSORTX_EMAIL = old_email)
  if (old_token != "") Sys.setenv(CIBERSORTX_TOKEN = old_token)
})

test_that("deconvolute errors on invalid algorithm", {
  bulk_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("Sample1", "Sample2")
    )
  )

  signature <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("CellType1", "CellType2")
    )
  )

  expect_error(
    deconvolute(
      algorithm = "invalid_algorithm",
      data = bulk_matrix,
      signature = signature
    ),
    "Unknown algorithm"
  )
})
