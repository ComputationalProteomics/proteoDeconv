#' Create signature matrix using CIBERSORTx
#'
#' Generates a signature matrix from reference profiles and cell type classes
#' using the CIBERSORTx Docker image. This signature matrix contains cell
#' type-specific marker proteins that will be used for deconvolution.
#'
#' @param refsample A numeric matrix containing reference profiles with genes as
#'   row names and samples as columns. This should be preprocessed data from
#'   pure cell populations.
#' @param phenoclasses A numeric matrix, data frame, or tibble containing cell
#'   type classification (0/1/2). If provided as a matrix or base data frame,
#'   the cell types are assumed to be in the row names. If provided as a tibble,
#'   the cell type identifiers should be included as a column. This is typically
#'   created using the create_phenoclasses() function.
#' @param g_min Minimum number of genes per cell type in the signature matrix.
#'   Default is 200.
#' @param g_max Maximum number of genes per cell type in the signature matrix.
#'   Default is 400.
#' @param q_value Q-value threshold for differential expression. Default is
#'   0.01.
#' @param replicates Number of replicates to use for building the reference file
#'   (only relevant when single_cell=TRUE). Default is 5.
#' @param sampling Fraction of available gene expression profiles selected by
#'   random sampling (only relevant when single_cell=TRUE). Default is 0.5.
#' @param fraction Fraction of cells of the same identity showing evidence of
#'   expression (only relevant when single_cell=TRUE). Default is 0.75.
#' @param filter Logical indicating whether to remove non-hematopoietic genes.
#'   Default is FALSE.
#' @param verbose Logical indicating whether to print detailed output. Default
#'   is FALSE.
#' @param QN Logical indicating whether to run quantile normalization. Default
#'   is FALSE.
#' @param single_cell Logical indicating whether to create signature from
#'   scRNA-Seq data. Default is FALSE.
#' @param use_sudo Logical indicating whether to use sudo for Docker commands.
#'   Default is FALSE.
#' @param ... Additional arguments passed to the function.
#'
#' @return A numeric matrix with genes as rows and cell types as columns,
#'   representing the expression profile of each cell type. This signature
#'   matrix can be used as input for deconvolution algorithms to estimate cell
#'   type proportions in mixed samples.
#'
#' @details This function uses the CIBERSORTx Docker image to construct a
#'   signature matrix based on the input reference profiles and cell type
#'   classifications. The CIBERSORTx token and email must be set as environment
#'   variables either in the project directory's .Renviron file or in the user's
#'   home directory .Renviron file.
#'
#' @examples
#' \dontrun{
#' # Load preprocessed pure samples data
#' pure_samples <- readRDS(system.file("extdata", "pure_samples_matrix.rds",
#'                                     package = "proteoDeconv"))
#'
#' # Create phenoclasses for the samples
#' mapping_rules <- list(
#'   "CD8+ T cells" = "CD8",
#'   "Monocytes" = "Mono"
#' )
#' phenoclasses <- create_phenoclasses(
#'   data = pure_samples,
#'   mapping_rules = mapping_rules
#' )
#'
#' # Create signature matrix
#' signature_matrix <- create_signature_matrix(
#'   refsample = pure_samples,
#'   phenoclasses = phenoclasses,
#'   g_min = 200,
#'   g_max = 400,
#'   q_value = 0.01
#' )
#' }
#'
#' @seealso \code{\link{create_phenoclasses}} for creating the phenoclasses
#'   input and \code{\link{deconvolute}} for using the signature matrix in
#'   deconvolution.
#'
#' @export
create_signature_matrix <- function(
  refsample = NULL,
  phenoclasses = NULL,
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
  use_sudo = FALSE,
  ...
) {
  if (!is.null(refsample) && !is.matrix(refsample)) {
    stop("refsample must be a matrix")
  }
  if (!is.null(phenoclasses)) {
    if (!(is.matrix(phenoclasses) || is.data.frame(phenoclasses))) {
      stop("phenoclasses must be a matrix or a tibble/data frame")
    }
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

  refsample_df <- NULL
  if (!is.null(refsample)) {
    refsample_df <- tibble::as_tibble(refsample, rownames = "Genes")
  }

  phenoclasses_df <- NULL
  if (!is.null(phenoclasses)) {
    if (tibble::is_tibble(phenoclasses)) {
      phenoclasses_df <- phenoclasses
    } else if (is.matrix(phenoclasses)) {
      phenoclasses_df <- as.data.frame(phenoclasses)
      phenoclasses_df <- cbind(rownames(phenoclasses_df), phenoclasses_df)
    } else if (is.data.frame(phenoclasses)) {
      if (
        is.null(rownames(phenoclasses)) ||
          all(rownames(phenoclasses) == as.character(1:nrow(phenoclasses)))
      ) {
        phenoclasses_df <- phenoclasses
      } else {
        phenoclasses_df <- cbind(
          rownames(phenoclasses),
          as.data.frame(phenoclasses)
        )
      }
    }
  }

  withr::with_tempdir({
    input_dir <- tempdir()
    output_dir <- tempdir()

    refsample_file <- if (!is.null(refsample_df)) {
      tempfile(tmpdir = input_dir, fileext = ".txt")
    } else {
      NULL
    }
    phenoclasses_file <- if (!is.null(phenoclasses_df)) {
      tempfile(tmpdir = input_dir, fileext = ".txt")
    } else {
      NULL
    }

    if (!is.null(refsample_df)) readr::write_tsv(refsample_df, refsample_file)
    if (!is.null(phenoclasses_df)) {
      utils::write.table(
        phenoclasses_df,
        phenoclasses_file,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
      )
    }

    docker_command <- glue::glue(
      "{if (use_sudo) 'sudo ' else ''}docker run --rm -v {input_dir}:/src/data:z -v {output_dir}:/src/outdir:z cibersortx/fractions ",
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
      "{ifelse(single_cell, '--single_cell TRUE', '')}"
    )

    if (verbose) {
      cat("Docker command:\n", docker_command, "\n")
    }

    command_output <- system(docker_command)
    if (command_output != 0) {
      stop(glue::glue("CIBERSORTx failed. Error code: {command_output}."))
    }

    refsample_id <- sub("\\.txt$", "", basename(refsample_file))
    signature_matrix_file <- grep(
      paste0("CIBERSORTx_", refsample_id),
      list.files(output_dir, pattern = "\\.txt$", full.names = TRUE),
      value = TRUE,
      fixed = TRUE
    )[1]

    if (is.na(signature_matrix_file)) {
      stop("No signature matrix file found in the output directory.")
    }

    signature_matrix_df <- readr::read_tsv(
      signature_matrix_file,
      show_col_types = FALSE
    )

    gene_col <- names(signature_matrix_df)[1]
    genes <- signature_matrix_df[[gene_col]]
    mat <- as.matrix(signature_matrix_df[, -1])
    rownames(mat) <- genes
    return(mat)
  })
}
