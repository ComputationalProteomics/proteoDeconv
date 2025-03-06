# Script to generate the signature matrix for the vignette
# This script should be run manually when updating the package data
# It requires Docker and CIBERSORTx credentials to be set up

library(proteoDeconv)
library(readr)
library(dplyr)


if (
  Sys.getenv("CIBERSORTX_EMAIL") == "" || Sys.getenv("CIBERSORTX_TOKEN") == ""
) {
  stop("Please set CIBERSORTX_EMAIL and CIBERSORTX_TOKEN environment variables")
}


pure_data_file <- system.file(
  "extdata",
  "pure_samples_matrix.rds",
  package = "proteoDeconv"
)
pure_samples <- readRDS(pure_data_file)


pure_samples <- pure_samples |>
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
  pure_samples,
  mapping_rules,
  verbose = TRUE
)

signature_matrix <- create_signature_matrix(
  refsample = pure_samples,
  phenoclasses = phenoclasses,
  g_min = 200,
  g_max = 400,
  q_value = 0.01,
  verbose = TRUE
)


saveRDS(pure_samples, "inst/extdata/pure_samples_processed.rds")
saveRDS(signature_matrix, "inst/extdata/cd8t_mono_signature_matrix.rds")
