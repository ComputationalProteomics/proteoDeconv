#' Deconvolute mixture data using CIBERSORT
#'
#' Applies the CIBERSORT algorithm to deconvolute bulk proteome data into constituent cell types
#' using a signature matrix. The function requires the original CIBERSORT.R script to be sourced.
#'
#' @param data A numeric matrix containing mixture data with genes as row names and samples as columns.
#' @param signature A numeric matrix containing signature data with genes as row names and cell types as columns.
#' @param QN Logical indicating whether quantile normalization is performed (default FALSE).
#' @param absolute Logical indicating whether an absolute score is computed (default FALSE).
#' @param abs_method Method for absolute scoring if absolute is TRUE (default "sig.score").
#' @param ... Additional arguments passed to the CIBERSORT function.
#'
#' @return A numeric matrix with samples as rows and cell types as columns, representing the estimated
#'         proportion of each cell type in each sample. The returned matrix excludes CIBERSORT's diagnostic
#'         columns (RMSE, P-value, Correlation).
#'
#' @details This function requires the original CIBERSORT.R script from the CIBERSORT website
#' (https://cibersortx.stanford.edu/) to be sourced before use. It writes temporary files for
#' the mixture data and signature matrix, calls the CIBERSORT function, and processes the results.
#'
#' @examples
#' \dontrun{
#' # First source the CIBERSORT.R script
#' source("/path/to/CIBERSORT.R")
#'
#' # Load example data and signature matrix
#' data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
#' mixed_samples <- readRDS(data_file)
#'
#' signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
#' signature_matrix <- readRDS(signature_file)
#'
#' # Run deconvolution
#' result <- deconvolute_cibersort(
#'   data = mixed_samples,
#'   signature = signature_matrix,
#'   QN = FALSE,
#'   absolute = FALSE
#' )
#' }
#'
#' @seealso \code{\link{deconvolute}} for a unified interface to multiple deconvolution methods.
#'
#' @export
deconvolute_cibersort <- function(
  data,
  signature,
  QN = FALSE,
  absolute = FALSE,
  abs_method = "sig.score",
  ...
) {
  if (!is.matrix(data)) {
    stop("Input 'data' must be a matrix")
  }
  if (!is.matrix(signature)) {
    stop("Input 'signature' must be a matrix")
  }

  if (!exists("CIBERSORT", mode = "function")) {
    stop(
      "Function 'CIBERSORT' not found. Please source the CIBERSORT.R script before running this function. You can download the script from the CIBERSORT website."
    )
  }

  cibersort_result <- withr::with_tempdir({
    expr_file <- tempfile()
    sig_file <- tempfile()

    expr_tbl <- tibble::as_tibble(data, rownames = "GeneSymbol")
    sig_tbl <- tibble::as_tibble(signature, rownames = "GeneSymbol")

    readr::write_tsv(expr_tbl, file = expr_file)
    readr::write_tsv(sig_tbl, file = sig_file)

    extras <- rlang::dots_list(
      sig_file,
      expr_file,
      perm = 0,
      QN = QN,
      absolute = absolute,
      abs_method = abs_method,
      ...,
      .homonyms = "last"
    )

    cibersort_call <- rlang::call2(CIBERSORT, !!!extras)
    output <- eval(cibersort_call)

    pruned_output <- output %>%
      .[, !colnames(.) %in% c("RMSE", "P-value", "Correlation")]

    return(pruned_output)
  })

  result <- as.matrix(cibersort_result)
  return(result)
}
