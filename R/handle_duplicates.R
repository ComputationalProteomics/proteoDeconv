#' Handle duplicate identifiers in an expression matrix
#'
#' Resolves duplicate row identifiers in an expression matrix using the specified method.
#'
#' @param data A numeric matrix containing expression data with identifiers as row names
#'        and samples as columns.
#' @param duplicate_mode A string specifying the approach to handle duplicates:
#'   \itemize{
#'     \item \code{"slice"}: Keep only the row with the maximum median value for each identifier (default)
#'     \item \code{"merge"}: Merge duplicate rows by taking the column-wise median of values
#'   }
#'
#' @return A numeric matrix with unique identifiers as row names. The number of rows will be
#'         equal to the number of unique identifiers in the input matrix.
#'
#'
#'
#' @examples
#' # Create example matrix with duplicate identifiers
#' mat <- matrix(1:12, nrow = 4, ncol = 3)
#' rownames(mat) <- c("ID1", "ID2", "ID1", "ID3")
#' colnames(mat) <- c("Sample1", "Sample2", "Sample3")
#'
#' # View original matrix
#' print(mat)
#'
#' # Handle duplicates by keeping rows with maximum median (default)
#' result1 <- handle_duplicates(mat)
#' print(result1)
#'
#' # Handle duplicates by merging rows
#' result2 <- handle_duplicates(mat, duplicate_mode = "merge")
#' print(result2)
#'
#' @export
handle_duplicates <- function(data, duplicate_mode = "slice") {
  if (!is.matrix(data)) {
    stop("Input must be a matrix")
  }
  if (is.null(rownames(data))) {
    stop("Matrix must have row names as identifiers")
  }

  ids <- rownames(data)
  missing_ids <- is.na(ids) | ids == ""
  if (any(missing_ids)) {
    ids[missing_ids] <- paste0(
      "UNKNOWN_",
      seq_along(ids)[missing_ids]
    )
    rownames(data) <- ids
  }

  if (!any(duplicated(ids))) {
    return(data)
  }

  unique_ids <- unique(ids)

  result <- matrix(NA, nrow = length(unique_ids), ncol = ncol(data))
  colnames(result) <- colnames(data)

  result <- switch(
    duplicate_mode,
    slice = {
      for (i in seq_along(unique_ids)) {
        id <- unique_ids[i]
        rows <- which(ids == id)
        if (length(rows) == 1) {
          result[i, ] <- data[rows, ]
        } else {
          medians <- apply(data[rows, , drop = FALSE], 1, function(row) {
            if (all(is.na(row))) -Inf else median(row, na.rm = TRUE)
          })
          best_row <- rows[which.max(medians)]
          result[i, ] <- data[best_row, ]
        }
      }
      result
    },
    merge = {
      for (i in seq_along(unique_ids)) {
        id <- unique_ids[i]
        rows <- which(ids == id)
        if (length(rows) == 1) {
          result[i, ] <- data[rows, ]
        } else {
          result[i, ] <- apply(
            data[rows, , drop = FALSE],
            2,
            median,
            na.rm = TRUE
          )
        }
      }
      result
    },
    stop("Unsupported duplicate mode: ", duplicate_mode)
  )

  rownames(result) <- unique_ids
  result
}
