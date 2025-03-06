#' Deconvolute bulk data using EPIC
#'
#' Runs the EPIC algorithm on preprocessed proteomic data to estimate cell type proportions
#' in mixed samples using a reference signature matrix.
#'
#' @param data A numeric matrix of the bulk proteome data with gene identifiers as row names
#'        and samples as columns.
#' @param signature A numeric matrix of reference signature profiles with gene identifiers as row names
#'        and cell types as columns.
#' @param with_other_cells Logical; if TRUE, EPIC will include an "other cells" component in the output
#'        to account for cell types not present in the reference. Default is TRUE.
#'
#' @return A numeric matrix with samples as rows and cell types as columns, representing the estimated
#'         proportion of each cell type in each sample.
#'
#' @details The function normalizes both the input data matrix and signature matrix using the
#'          `handle_scaling` function before running the EPIC deconvolution. The EPIC algorithm
#'          uses constrained least squares regression to estimate cell type proportions.
#'
#' @examples
#' \dontrun{
#' # Load example data and signature matrix
#' data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
#' mixed_samples <- readRDS(data_file)
#'
#' signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
#' signature_matrix <- readRDS(signature_file)
#'
#' # Run EPIC deconvolution
#' result <- deconvolute_epic(
#'   data = mixed_samples,
#'   signature = signature_matrix,
#'   with_other_cells = TRUE
#' )
#'
#' # View first few rows of the result
#' head(result)
#' }
#'
#' @seealso \code{\link{handle_scaling}} for preprocessing inputs before deconvolution,
#'          \code{\link{deconvolute}} for a unified interface to multiple deconvolution methods.
#'
#' @references Racle, J. et al. (2017). Simultaneous enumeration of cancer and immune cell types
#'             from bulk tumor gene expression data. eLife, 6:e26476.
#'
#' @export
deconvolute_epic <- function(
  data,
  signature,
  with_other_cells = TRUE
) {
  if (!is.matrix(data)) {
    stop("data must be a matrix")
  }
  if (!is.matrix(signature)) {
    stop("signature must be a matrix")
  }
  if (is.null(rownames(data)) || is.null(rownames(signature))) {
    stop("Both matrices must have gene identifiers as row names")
  }

  signature_mat_epic <- handle_scaling(
    signature,
    tpm = TRUE,
    unlog = FALSE
  )

  preprocessed_mat_epic <- handle_scaling(
    data,
    tpm = TRUE,
    unlog = FALSE
  )

  epic_signature <- list(
    refProfiles = signature_mat_epic,
    sigGenes = rownames(signature_mat_epic)
  )

  epic_res <- suppressWarnings(EPIC::EPIC(
    preprocessed_mat_epic,
    withOtherCells = with_other_cells,
    reference = epic_signature
  ))

  return(epic_res$mRNAProportions)
}
