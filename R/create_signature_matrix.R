#' Create Signature Matrix using CIBERSORTx
#'
#' This function generates a signature matrix from a reference sample and cell type classes using the CIBERSORTx Docker tool.
#'
#' @param refsample A data frame containing the reference profiles (with replicates). Required.
#' @param phenoclasses A data frame containing the cell type classes. Required.
#' @param g_min Minimum number of genes per cell type in the signature matrix. Default is 200.
#' @param g_max Maximum number of genes per cell type in the signature matrix. Default is 400.
#' @param q_value Q-value threshold for differential expression. Default is 0.01.
#' @param replicates Number of replicates to use for building the reference file. Default is 5.
#' @param sampling Fraction of available gene expression profiles selected by random sampling. Default is 0.5.
#' @param fraction Fraction of cells of the same identity showing evidence of expression. Default is 0.75.
#' @param filter Logical indicating whether to remove non-hematopoietic genes. Default is FALSE.
#' @param verbose Logical indicating whether to print detailed output. Default is FALSE.
#' @param QN Logical indicating whether to run quantile normalization. Default is FALSE.
#' @param single_cell Logical indicating whether to create signature from scRNA-Seq data. Default is FALSE.
#' @param ... Additional arguments passed to the function.
#'
#' @return A data frame containing the generated signature matrix.
#'
#' @details This function uses the CIBERSORTx Docker image to construct a signature matrix based on the input reference profiles and cell type classifications.
#'
#' @importFrom readr write_tsv read_tsv
#' @importFrom glue glue
#' @importFrom withr with_tempdir
#' @export
create_signature_matrix <- function(
    refsample = NULL, phenoclasses = NULL,
    g_min = 200,
    g_max = 400,
    q_value = 0.01,
    replicates = 5,
    sampling = 0.5,
    fraction = 0.75,
    filter = FALSE,
    verbose = FALSE,
    QN = FALSE,
    single_cell = FALSE,
    ...) {
  username <- Sys.getenv("CIBERSORTX_EMAIL")
  token <- Sys.getenv("CIBERSORTX_TOKEN")
  if (username == "" || token == "") {
    stop("CIBERSORTX_EMAIL and CIBERSORTX_TOKEN environment variables must be set.")
  }
  if (system("docker --version", intern = TRUE, ignore.stderr = TRUE) == 127) {
    stop("Docker is not installed or not running.")
  }

  withr::with_tempdir({
    input_dir <- tempdir()
    output_dir <- tempdir()

    refsample_file <- if (!is.null(refsample)) tempfile(tmpdir = input_dir) else NULL
    phenoclasses_file <- if (!is.null(phenoclasses)) tempfile(tmpdir = input_dir) else NULL

    if (!is.null(refsample)) readr::write_tsv(refsample, refsample_file)
    if (!is.null(phenoclasses)) readr::write_tsv(phenoclasses, phenoclasses_file, col_names = FALSE)

    docker_command <- glue::glue(
      "{if (Sys.info()[['sysname']] == 'Linux') 'sudo ' else ''}docker run --rm -v {input_dir}:/src/data:z -v {output_dir}:/src/outdir:z cibersortx/fractions ",
      "--verbose {ifelse(verbose, 'TRUE', 'FALSE')} ",
      "--username {username} ",
      "--token {token} ",
      "{if (!is.null(refsample_file)) paste('--refsample', refsample_file) else ''} ",
      "{if (!is.null(phenoclasses_file)) paste('--phenoclasses', phenoclasses_file) else ''} ",
      "--G.min {g_min} ",
      "--G.max {g_max} ",
      "--q.value {q_value} ",
      "--filter {ifelse(filter, 'TRUE', 'FALSE')} ",
      "--QN {ifelse(QN, 'TRUE', 'FALSE')} ",
      "--k.max 999 ",
      "--remake FALSE ",
      "--replicates {replicates} ",
      "--sampling {sampling} ",
      "--fraction {fraction} ",
      "--single_cell {ifelse(single_cell, 'TRUE', 'FALSE')}"
    )
    if (verbose) {
      cat("Docker command:\n", docker_command, "\n")
    }
    command_output <- system(docker_command)
    if (command_output != 0) {
      stop(glue::glue("CIBERSORTx failed. Error code: {command_output}."))
    }
    signature_matrix_file <- grep(basename(refsample_file), list.files(output_dir, pattern = "\\.txt$", full.names = TRUE), value = TRUE)
    if (length(signature_matrix_file) == 0) {
      stop("No signature matrix file found in the output directory.")
    }
    signature_matrix <- readr::read_tsv(signature_matrix_file, show_col_types = FALSE)
    signature_matrix
  })
}
