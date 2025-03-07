#' Handle missing values
#'
#' Imputes missing values in an expression matrix using various methods,
#' including simple replacement with the lowest non-zero value or more
#' sophisticated imputation techniques via the MsCoreUtils package.
#'
#' @param data A numeric matrix containing the data to be processed, with
#'   identifiers as row names and samples as columns.
#' @param imputation_mode A string specifying the imputation method to use:
#'   \itemize{
#'     \item \code{"lowest_value"} (default): Replaces NAs with the lowest non-zero value in the matrix
#'     \item \code{"knn"}: k-nearest neighbor imputation
#'     \item \code{"zero"}: Replace NAs with zeros
#'     \item \code{"MLE"}: Maximum likelihood estimation
#'     \item \code{"bpca"}: Bayesian principal component analysis
#'     \item \code{"RF"}: Random Forest imputation
#'     \item \code{"min"}: Replace NAs with minimum value in each column
#'     \item \code{"MinDet"}: Deterministic minimum value imputation
#'     \item \code{"MinProb"}: Probabilistic minimum value imputation
#'     \item \code{"QRILC"}: Quantile regression imputation of left-censored data
#'     \item \code{"mixed"}: Mixed imputation based on feature-wise missingness
#'     \item \code{"nbavg"}: Impute with average of neighbors
#'   }
#' @param ... Additional arguments passed to \code{MsCoreUtils::impute_matrix}.
#'   See the documentation of that function for method-specific parameters.
#'
#' @return A matrix with the same dimensions as the input, but with missing
#'   values imputed according to the specified method.
#'
#'
#' @examples
#' # Create example matrix with missing values
#' mat <- matrix(c(1.2, 3.4, NA, 5.6, NA, 7.8, 9.0, 2.1, 4.3), nrow = 3, ncol = 3)
#' rownames(mat) <- c("Protein1", "Protein2", "Protein3")
#' colnames(mat) <- c("Sample1", "Sample2", "Sample3")
#'
#' # View original matrix
#' print(mat)
#'
#' # Impute missing values with lowest non-zero value (default)
#' result1 <- handle_missing_values(mat)
#' print(result1)
#'
#' # Impute missing values using k-nearest neighbors
#' \dontrun{
#' # Requires the 'impute' package
#' result2 <- handle_missing_values(mat, imputation_mode = "knn")
#' print(result2)
#' }
#'
#' @seealso \code{\link[MsCoreUtils]{impute_matrix}} for detailed description of
#'   the imputation methods
#' @export

handle_missing_values <- function(
  data,
  imputation_mode = "lowest_value",
  ...
) {
  if (!is.matrix(data)) {
    stop("Input must be a matrix")
  }

  row_names <- rownames(data)
  col_names <- colnames(data)

  result <- switch(
    imputation_mode,
    "lowest_value" = {
      lowest_value <- min_nonzero(data)
      data[is.na(data)] <- lowest_value
      data
    },
    "knn" = MsCoreUtils::impute_matrix(data, method = "knn", ...),
    "zero" = MsCoreUtils::impute_matrix(data, method = "zero", ...),
    "MLE" = MsCoreUtils::impute_matrix(data, method = "MLE", ...),
    "bpca" = MsCoreUtils::impute_matrix(data, method = "bpca", ...),
    "RF" = MsCoreUtils::impute_matrix(data, method = "RF", ...),
    "min" = MsCoreUtils::impute_matrix(data, method = "min", ...),
    "MinDet" = MsCoreUtils::impute_matrix(data, method = "MinDet", ...),
    "MinProb" = MsCoreUtils::impute_matrix(data, method = "MinProb", ...),
    "QRILC" = MsCoreUtils::impute_matrix(data, method = "QRILC", ...),
    "mixed" = MsCoreUtils::impute_matrix(data, method = "mixed", ...),
    "nbavg" = MsCoreUtils::impute_matrix(data, method = "nbavg", ...),
    stop(paste("Unsupported imputation mode:", imputation_mode))
  )

  rownames(result) <- row_names
  colnames(result) <- col_names

  return(result)
}

min_nonzero <- function(x) {
  non_zero_values <- x[x != 0]
  min(non_zero_values, na.rm = TRUE)
}
