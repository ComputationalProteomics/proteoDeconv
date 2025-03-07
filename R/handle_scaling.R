#' Handle scaling of protein abundance data
#'
#' Scales protein abundance data by optionally converting from log2 scale to
#' linear scale and/or applying TPM-like normalization for proteomics data
#' comparison.
#'
#' @param data A numeric matrix containing protein abundance data with
#'   identifiers as row names and samples as columns.
#' @param unlog A logical value indicating whether to unlog the data (convert
#'   from log2 scale). Default is FALSE.
#' @param tpm A logical value indicating whether to apply TPM-like normalization
#'   (adapting Transcripts Per Million for proteomics). Default is FALSE.
#'
#' @return A numeric matrix with the scaled protein abundance data according to
#'   the specified options.
#'
#' @details This function combines the functionality of
#'   \code{\link{unlog2_data}} and \code{\link{convert_to_tpm}} to provide a
#'   flexible way to handle common scaling operations for proteomics data.
#'
#' @examples
#' # Create log2-transformed protein abundance matrix
#' log2_mat <- matrix(rnorm(12, mean = 10, sd = 2), nrow = 4, ncol = 3)
#' rownames(log2_mat) <- paste0("Protein", 1:4)
#' colnames(log2_mat) <- paste0("Sample", 1:3)
#'
#' # Convert from log2 to linear scale only
#' linear_mat <- handle_scaling(log2_mat, unlog = TRUE, tpm = FALSE)
#'
#' # Convert from log2 to linear scale and then apply TPM-like normalization
#' tpm_mat <- handle_scaling(log2_mat, unlog = TRUE, tpm = TRUE)
#'
#' # Apply TPM-like normalization to already linear protein abundance data
#' linear_data <- matrix(abs(rnorm(12, mean = 500, sd = 200)), nrow = 4, ncol = 3)
#' tpm_only <- handle_scaling(linear_data, unlog = FALSE, tpm = TRUE)
#'
#' @seealso \code{\link{unlog2_data}} for just converting from log2 scale,
#'   \code{\link{convert_to_tpm}} for just applying TPM-like normalization
#'
#' @export
handle_scaling <- function(
  data,
  unlog = FALSE,
  tpm = FALSE
) {
  if (!is.matrix(data)) {
    stop("Input must be a matrix")
  }

  row_names <- rownames(data)
  col_names <- colnames(data)

  processed_mat <- data

  if (unlog) {
    processed_mat <- unlog2_data(processed_mat)
  }

  if (tpm) {
    processed_mat <- convert_to_tpm(processed_mat)
  }

  return(processed_mat)
}

#' Convert data from log2 scale to linear scale
#'
#' Transforms log2-transformed expression data back to linear scale by calculating 2^x
#' for each value in the input matrix.
#'
#' @param data A numeric matrix containing log2-transformed data with identifiers as row names
#'        and samples as columns.
#'
#' @return A numeric matrix with the same dimensions as the input, with values converted to
#'         linear scale (2^x).
#'
#' @details This function can be used to convert log2-transformed proteomics or gene expression
#' data back to their original scale for downstream analyses that require linear values.
#' It preserves the row and column names of the input matrix.
#'
#' @examples
#' # Create log2-transformed data matrix
#' log2_mat <- matrix(rnorm(12, mean = 10, sd = 2), nrow = 4, ncol = 3)
#' rownames(log2_mat) <- paste0("Protein", 1:4)
#' colnames(log2_mat) <- paste0("Sample", 1:3)
#'
#' # View original log2 values
#' print(log2_mat)
#'
#' # Convert to linear scale
#' linear_mat <- unlog2_data(log2_mat)
#' print(linear_mat)
#'
#' @export
unlog2_data <- function(data) {
  if (!is.matrix(data)) {
    stop("Input must be a matrix")
  }

  row_names <- rownames(data)
  col_names <- colnames(data)

  processed_mat <- 2^data

  rownames(processed_mat) <- row_names
  colnames(processed_mat) <- col_names

  return(processed_mat)
}

#' Convert protein abundance data to TPM-like normalization
#'
#' Normalizes protein abundance data using a TPM-like approach (Transcripts Per Million),
#' adapting this RNA-seq normalization method for use with proteomics data.
#'
#' @param data A numeric matrix containing protein abundance data with identifiers as row names
#'        and samples as columns. Values should be in linear scale, not log-transformed.
#'
#' @return A numeric matrix with the same dimensions as the input, with values normalized
#'         to a TPM-like scale (sum of each column equals 1 million).
#'
#' @details This function applies a TPM-like normalization to proteomics data, where each protein
#' abundance value is scaled by the total abundance in the sample and multiplied by 1 million.
#' While TPM was originally designed for RNA-seq data (where it accounts for transcript length),
#' this adapted version provides a similar relative scaling benefit for proteomics data by
#' normalizing for differences in total protein abundance between samples.
#'
#' This normalization enables more meaningful comparisons of protein abundance levels
#' across samples with different total detected protein amounts. The resulting values
#' represent relative abundances where the sum for each sample equals 1 million.
#'
#' If your input data is log-transformed, use \code{\link{unlog2_data}} first to convert it to
#' linear scale before applying this normalization.
#'
#' @examples
#' # Create example protein abundance data matrix
#' prot_mat <- matrix(abs(rnorm(12, mean = 500, sd = 200)), nrow = 4, ncol = 3)
#' rownames(prot_mat) <- paste0("Protein", 1:4)
#' colnames(prot_mat) <- paste0("Sample", 1:3)
#'
#' # View original values and column sums
#' print(prot_mat)
#' print(colSums(prot_mat))
#'
#' # Convert to TPM-like normalization
#' tpm_mat <- convert_to_tpm(prot_mat)
#'
#' # Verify that column sums equal 1 million
#' print(tpm_mat)
#' print(colSums(tpm_mat))
#'
#' @export
convert_to_tpm <- function(data) {
  if (!is.matrix(data)) {
    stop("Input must be a matrix")
  }

  if (any(data < 0, na.rm = TRUE)) {
    warning(
      "Negative values detected. Input data should be in linear scale, not log-transformed."
    )
  } else if (max(data, na.rm = TRUE) < 100) {
    warning(
      "Maximum value is small. Input data might be log-transformed. Consider using unlog2_data() first."
    )
  }

  row_names <- rownames(data)
  col_names <- colnames(data)

  processed_mat <- data
  col_sums <- colSums(processed_mat, na.rm = TRUE)

  for (i in seq_along(col_sums)) {
    if (col_sums[i] > 0) {
      processed_mat[, i] <- (processed_mat[, i] / col_sums[i]) * 1e6
    }
  }

  rownames(processed_mat) <- row_names
  colnames(processed_mat) <- col_names

  return(processed_mat)
}
