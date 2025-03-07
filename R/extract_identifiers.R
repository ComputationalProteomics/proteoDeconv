#' Extract primary identifiers from e.g. protein groups
#'
#' Extracts the primary identifier from compound row identifiers in a matrix by
#' splitting on a separator and keeping only the first entry.
#'
#' @param data A numeric matrix with compound identifiers (e.g. protein/gene
#'   groups) as row names.
#' @param separator A character string used to separate multiple identifiers.
#'   Default is ";".
#'
#' @return A matrix with simplified identifiers as row names, where each
#'   identifier is the first element from the original compound identifier.
#'
#'
#' @examples
#' # Create matrix with compound identifiers (like protein groups)
#' mat <- matrix(1:9, nrow = 3, ncol = 3)
#' rownames(mat) <- c("P04637;P02340", "Q15796;O35182", "P01308;P01315;P01317")
#' colnames(mat) <- c("Sample1", "Sample2", "Sample3")
#'
#' # View original matrix
#' print(mat)
#'
#' # Extract primary identifiers
#' result <- extract_identifiers(mat)
#' print(result)
#'
#' @export
extract_identifiers <- function(data, separator = ";") {
  if (!is.matrix(data)) {
    stop("Input must be a matrix")
  }

  if (is.null(rownames(data))) {
    stop("Matrix must have row names as identifiers")
  }

  ids <- rownames(data)

  simplified_ids <- strsplit(ids, separator) |>
    lapply(function(x) x[1]) |>
    unlist()

  rownames(data) <- simplified_ids

  return(data)
}
