#' Create Signature Matrix using CIBERSORTx
#'
#' This function generates a signature matrix from a reference sample and phenoclasses using the CIBERSORTx tool.
#'
#' @param refsample A data frame containing the reference profiles (with replicates). This parameter is required.
#' @param phenoclasses A data frame containing the cell type classes. This parameter is required.
#' @param g_min Minimum number of genes per cell type in the signature matrix. Default is 50.
#' @param g_max Maximum number of genes per cell type in the signature matrix. Default is 150.
#' @param q_value Q-value threshold for differential expression. Default is 0.3.
#' @param replicates Number of replicates to use for building the reference file. Default is 5.
#' @param sampling Fraction of available gene expression profiles selected using random sampling. Default is 0.5.
#' @param fraction Fraction of cells of the same identity showing evidence of expression. Default is 0.75.
#' @param filter Logical indicating whether to remove non-hematopoietic genes. Default is FALSE.
#' @param verbose Logical indicating whether to print detailed output. Default is TRUE.
#' @param ... Additional arguments passed to the function.
#'
#' @return A data frame containing the signature matrix.
#'
#' @details This function uses the CIBERSORTx Docker image to generate a signature matrix.
#'
#' @importFrom readr write_tsv read_tsv
#' @importFrom glue glue
#' @importFrom withr with_tempdir
#' @export
create_signature_matrix <- function(
    refsample, phenoclasses,
    g_min = 50,
    g_max = 150,
    q_value = 0.3,
    replicates = 5,
    sampling = 0.5,
    fraction = 0.75,
    filter = FALSE,
    verbose = TRUE, ...) {
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

    refsample_file <- tempfile(tmpdir = input_dir)
    phenoclasses_file <- tempfile(tmpdir = input_dir)

    readr::write_tsv(refsample, refsample_file)
    readr::write_tsv(phenoclasses, phenoclasses_file)

    docker_command <- glue::glue(
      "docker run -v {input_dir}:/src/data:z -v {output_dir}:/src/outdir:z cibersortx/fractions ",
      "--verbose {ifelse(verbose, 'TRUE', 'FALSE')} ",
      "--username {username} ",
      "--token {token} ",
      "--refsample {refsample_file} ",
      "--phenoclasses {phenoclasses_file} ",
      "--G.min {g_min} ",
      "--G.max {g_max} ",
      "--q.value {q_value} ",
      "--filter {ifelse(filter, 'TRUE', 'FALSE')} ",
      "--k.max 999 ",
      "--remake FALSE ",
      "--replicates {replicates} ",
      "--sampling {sampling} ",
      "--fraction {fraction}"
    )
    if (verbose) {
      cat("Docker command:\n", docker_command, "\n")
    }

    command_output <- system(docker_command)

    if (command_output != 0) {
      stop(glue::glue("CIBERSORTx failed. Error code: {command_output}"))
    }

    signature_matrix_file <- grep(basename(refsample_file), list.files(output_dir, pattern = "\\.txt$", full.names = TRUE), value = TRUE)
    if (length(signature_matrix_file) == 0) {
      stop("No signature matrix file found in the output directory.")
    }

    signature_matrix <- readr::read_tsv(signature_matrix_file, show_col_types = FALSE)
    signature_matrix
  })
}
