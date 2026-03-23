skip_if_no_cibersortx_live <- function(message) {
  docker_is_available <- getFromNamespace("docker_is_available", "proteoDeconv")

  skip_if_not(
    identical(Sys.getenv("RUN_CIBERSORTX_LIVE"), "true") &&
      identical(Sys.info()[["sysname"]], "Linux") &&
      Sys.getenv("CIBERSORTX_EMAIL") != "" &&
      Sys.getenv("CIBERSORTX_TOKEN") != "" &&
      docker_is_available(),
    message
  )
}

read_test_matrix <- function(path) {
  skip_if(path == "", "Test data not available")

  data <- readRDS(path)
  if (is.data.frame(data)) {
    gene_col <- which(vapply(data, is.character, logical(1)))[1]
    genes <- data[[gene_col]]
    data <- as.matrix(data[, -gene_col])
    rownames(data) <- genes
  }

  data
}

load_test_pure_samples <- function() {
  read_test_matrix(system.file(
    "extdata",
    "pure_samples_matrix.rds",
    package = "proteoDeconv"
  ))
}

load_test_mixed_samples <- function() {
  read_test_matrix(system.file(
    "extdata",
    "mixed_samples_matrix.rds",
    package = "proteoDeconv"
  ))
}

load_test_signature_matrix <- function() {
  read_test_matrix(system.file(
    "extdata",
    "cd8t_mono_signature_matrix.rds",
    package = "proteoDeconv"
  ))
}

cibersortx_test_mapping_rules <- function() {
  list(
    "CD8+ T cells" = "CD8",
    "Monocytes" = "Mono"
  )
}
