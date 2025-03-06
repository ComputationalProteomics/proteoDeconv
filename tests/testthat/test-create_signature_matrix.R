test_that("create_signature_matrix validates inputs correctly", {
  expect_error(
    create_signature_matrix(
      refsample = data.frame(Genes = c("Gene1"), Sample1 = c(1)),
      phenoclasses = matrix(1, nrow = 1, ncol = 1)
    ),
    "refsample must be a matrix"
  )
})

test_that("create_signature_matrix checks for credentials", {
  original_email <- Sys.getenv("CIBERSORTX_EMAIL")
  original_token <- Sys.getenv("CIBERSORTX_TOKEN")

  Sys.unsetenv("CIBERSORTX_EMAIL")
  Sys.unsetenv("CIBERSORTX_TOKEN")

  expect_error(
    create_signature_matrix(
      refsample = matrix(
        1,
        nrow = 1,
        ncol = 1,
        dimnames = list("Gene1", "Sample1")
      ),
      phenoclasses = matrix(
        1,
        nrow = 1,
        ncol = 1,
        dimnames = list("CellType", "Sample1")
      )
    ),
    "CIBERSORTX_EMAIL and CIBERSORTX_TOKEN environment variables must be set"
  )

  if (original_email != "") Sys.setenv(CIBERSORTX_EMAIL = original_email)
  if (original_token != "") Sys.setenv(CIBERSORTX_TOKEN = original_token)
})

test_that("create_signature_matrix works with real data", {
  skip_if_not(
    Sys.getenv("CIBERSORTX_EMAIL") != "" &&
      Sys.getenv("CIBERSORTX_TOKEN") != "" &&
      system("docker --version", ignore.stderr = TRUE) == 0,
    "Skipping test - Docker or CIBERSORTx credentials not available"
  )

  pure_data_file <- system.file(
    "extdata",
    "pure_samples_matrix.rds",
    package = "proteoDeconv"
  )
  skip_if(pure_data_file == "", "Test data not available")

  pure_samples <- readRDS(pure_data_file)
  if (is.data.frame(pure_samples)) {
    gene_col <- which(sapply(pure_samples, is.character))[1]
    genes <- pure_samples[[gene_col]]
    mat_data <- as.matrix(pure_samples[, -gene_col])
    rownames(mat_data) <- genes
    pure_samples <- mat_data
  }

  processed_samples <- pure_samples |>
    extract_identifiers() |>
    handle_missing_values(imputation_mode = "lowest_value") |>
    handle_duplicates(duplicate_mode = "slice") |>
    handle_scaling(unlog = FALSE, tpm = TRUE)

  mapping_rules <- list(
    "CD8+ T cells" = "CD8",
    "Monocytes" = "Mono"
  )

  phenoclasses <- create_phenoclasses(
    data = processed_samples,
    mapping_rules = mapping_rules,
    verbose = TRUE
  )

  expected_cell_types <- names(mapping_rules)

  result <- create_signature_matrix(
    refsample = processed_samples,
    phenoclasses = phenoclasses,
    g_min = 5,
    g_max = 10,
    q_value = 0.3,
    replicates = 1,
    sampling = 1.0,
    fraction = 0.5,
    verbose = TRUE
  )

  expect_true(is.matrix(result))
  expect_true(nrow(result) > 0)
  expect_true(ncol(result) > 0)

  expect_true(any(expected_cell_types %in% colnames(result)))
})

test_that("create_signature_matrix matches reference workflow", {
  skip_if_not(
    Sys.getenv("CIBERSORTX_EMAIL") != "" &&
      Sys.getenv("CIBERSORTX_TOKEN") != "" &&
      system("docker --version", ignore.stderr = TRUE) == 0,
    "Skipping full workflow test - prerequisites not met"
  )

  pure_data_file <- system.file(
    "extdata",
    "pure_samples_matrix.rds",
    package = "proteoDeconv"
  )
  skip_if(pure_data_file == "", "Test data not available")

  pure_samples <- readRDS(pure_data_file)
  if (is.data.frame(pure_samples)) {
    gene_col <- which(sapply(pure_samples, is.character))[1]
    genes <- pure_samples[[gene_col]]
    mat_data <- as.matrix(pure_samples[, -gene_col])
    rownames(mat_data) <- genes
    pure_samples <- mat_data
  }

  processed_samples <- pure_samples |>
    extract_identifiers() |>
    update_gene_symbols(verbose = FALSE) |>
    handle_missing_values(imputation_mode = "lowest_value") |>
    handle_duplicates(duplicate_mode = "slice") |>
    handle_scaling(unlog = FALSE, tpm = TRUE)

  mapping_rules <- list(
    "CD8+ T cells" = "CD8",
    "Monocytes" = "Mono"
  )

  phenoclasses <- create_phenoclasses(
    data = processed_samples,
    mapping_rules = mapping_rules,
    verbose = TRUE
  )

  signature_matrix <- create_signature_matrix(
    refsample = processed_samples,
    phenoclasses = phenoclasses,
    g_min = 10,
    g_max = 20,
    q_value = 0.2,
    replicates = 1,
    verbose = TRUE
  )

  expect_true(is.matrix(signature_matrix))

  expected_columns <- names(mapping_rules)
  cell_types_in_result <- colnames(signature_matrix)

  expect_true(any(expected_columns %in% cell_types_in_result))

  expect_true(all(rownames(signature_matrix) %in% rownames(processed_samples)))
})
