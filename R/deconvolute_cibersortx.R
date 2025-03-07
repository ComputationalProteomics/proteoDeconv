#' Deconvolute using CIBERSORTx Docker image
#'
#' Performs deconvolution of bulk proteome data into constituent cell types
#' using the CIBERSORTx Docker image. This function handles the interaction with
#' the Docker container and processes the results.
#'
#' @param data A numeric matrix containing mixture data with genes as row names
#'   and samples as columns.
#' @param signature A numeric matrix containing the signature matrix with genes
#'   as row names and cell types as columns.
#' @param perm An integer specifying the number of permutations to be performed.
#'   Default is 1.
#' @param rmbatch_S_mode A logical value indicating whether to remove batch
#'   effects in source GEPs mode. Default is FALSE.
#' @param source_GEPs A matrix containing the source gene expression profiles.
#'   Required if \code{rmbatch_S_mode} is TRUE.
#' @param use_cibersortx A logical value indicating whether to use CIBERSORTx.
#'   Default is TRUE.
#' @param rmbatch_B_mode A logical value indicating whether to remove batch
#'   effects in bulk mode. Default is FALSE.
#' @param QN A logical value indicating whether to perform quantile
#'   normalization. Default is FALSE.
#' @param absolute A logical value indicating whether to use absolute mode.
#'   Default is FALSE.
#' @param abs_method A character string specifying the method to use for
#'   absolute mode. Default is "sig.score".
#' @param use_sudo A logical value indicating whether to use sudo for Docker
#'   commands. Default is FALSE.
#'
#' @return A numeric matrix with samples as rows and cell types as columns,
#'   representing the estimated proportion of each cell type in each sample.
#'
#' @details This function requires the CIBERSORTx Docker image to be installed
#' and the \code{CIBERSORTX_EMAIL} and \code{CIBERSORTX_TOKEN} environment
#' variables to be set. You can get these credentials by registering at the
#' CIBERSORTx website (https://cibersortx.stanford.edu/).
#'
#' The function creates temporary files for the mixture data and signature
#' matrix, runs the CIBERSORTx Docker container, and processes the results. Note
#' that absolute mode is not currently supported in the Docker version.
#'
#' @examples
#' \dontrun{
#' # Set required environment variables (ideally in .Renviron)
#' Sys.setenv(CIBERSORTX_EMAIL = "your.email@example.com")
#' Sys.setenv(CIBERSORTX_TOKEN = "your-token-here")
#'
#' # Load example data and signature matrix
#' data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
#' mixed_samples <- readRDS(data_file)
#'
#' signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
#' signature_matrix <- readRDS(signature_file)
#'
#' # Run deconvolution with CIBERSORTx Docker
#' result <- deconvolute_cibersortx(
#'   data = mixed_samples,
#'   signature = signature_matrix,
#'   perm = 100,
#'   QN = TRUE
#' )
#' }
#'
#' @seealso \code{\link{deconvolute_cibersort}} for using the R implementation
#'   of CIBERSORT, \code{\link{deconvolute}} for a unified interface to multiple
#'   deconvolution methods.
#'
#' @export
deconvolute_cibersortx <- function(
  data,
  signature,
  perm = 1,
  rmbatch_S_mode = FALSE,
  source_GEPs = NULL,
  use_cibersortx = TRUE,
  rmbatch_B_mode = FALSE,
  QN = FALSE,
  absolute = FALSE,
  abs_method = "sig.score",
  use_sudo = FALSE
) {
  if (absolute == TRUE) {
    stop(
      "Absolute mode is not supported in the CIBERSORTx Docker image. Use the online version instead, or use the regular CIBERSORT method."
    )
  }

  username <- Sys.getenv("CIBERSORTX_EMAIL")
  token <- Sys.getenv("CIBERSORTX_TOKEN")
  if (username == "" || token == "") {
    stop(
      "CIBERSORTX_EMAIL and CIBERSORTX_TOKEN environment variables must be set."
    )
  }

  if (system("docker --version", intern = TRUE, ignore.stderr = TRUE) == 127) {
    stop("Docker is not installed or not running.")
  }

  data_tibble <- tibble::as_tibble(data, rownames = "GeneSymbol")
  sig_tibble <- tibble::as_tibble(signature, rownames = "GeneSymbol")

  source_GEPs_tibble <- NULL
  if (rmbatch_S_mode) {
    if (is.null(source_GEPs)) {
      stop("source_GEPs is required when rmbatch_S_mode is TRUE")
    }
    if (!is.matrix(source_GEPs)) {
      stop("source_GEPs must be a matrix")
    }
    source_GEPs_tibble <- tibble::as_tibble(
      source_GEPs,
      rownames = "GeneSymbol"
    )
  }

  withr::with_tempdir({
    input_dir <- tempdir()
    output_dir <- tempdir()

    input_data_file <- tempfile(tmpdir = input_dir)
    signature_file <- tempfile(tmpdir = input_dir)

    source_GEPs_file <- NULL
    if (rmbatch_S_mode) {
      source_GEPs_file <- tempfile(tmpdir = input_dir)
      readr::write_tsv(source_GEPs_tibble, source_GEPs_file)
    }

    readr::write_tsv(data_tibble, input_data_file)
    readr::write_tsv(sig_tibble, signature_file)

    label <- uuid::UUIDgenerate(TRUE)

    docker_command <- glue::glue(
      "{if (use_sudo) 'sudo ' else ''}docker run --rm -v {input_dir}:/src/data:z -v {output_dir}:/src/outdir:z cibersortx/fractions ",
      "--verbose TRUE ",
      "--username {username} ",
      "--token {token} ",
      "--mixture {input_data_file} ",
      "--sigmatrix {signature_file} ",
      "--perm {perm} ",
      "--label {label} ",
      "--rmbatchBmode {ifelse(rmbatch_B_mode, 'TRUE', 'FALSE')} ",
      "--rmbatchSmode {ifelse(rmbatch_S_mode, 'TRUE', 'FALSE')} ",
      "--sourceGEPs {source_GEPs_file %||% signature_file} ",
      "--QN {ifelse(QN, 'TRUE', 'FALSE')} ",
      "--absolute {ifelse(absolute, 'TRUE', 'FALSE')} ",
      "--abs_method {abs_method} "
    )
    message("Docker command:\n", docker_command, "\n")

    command_output <- system(docker_command)

    if (command_output != 0) {
      stop(glue::glue("CIBERSORTx failed. Error code: {command_output}"))
    }

    result_file <- glue::glue("{output_dir}/CIBERSORTx_{label}_Results.txt")
    if (!file.exists(result_file)) {
      stop("No result file found in the output directory.")
    }

    cibersortx_result <- readr::read_tsv(result_file, show_col_types = FALSE)

    result_matrix <- convert_cibersortx_output_to_matrix(cibersortx_result)
    return(result_matrix)
  })
}

convert_cibersortx_output_to_matrix <- function(data) {
  extra_cols <- c("P.value", "Correlation", "RMSE", "P-value")

  if (!("Mixture" %in% colnames(data))) {
    stop("No Mixture column found")
  }

  cell_cols <- setdiff(colnames(data), c("Mixture", extra_cols))

  result_matrix <- as.matrix(data[, cell_cols])
  rownames(result_matrix) <- data$Mixture

  return(result_matrix)
}
