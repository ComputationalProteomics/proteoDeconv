#' Create phenoclasses matrix for cell type deconvolution
#'
#' Generates a phenotype classification matrix from sample data for use in cell
#' type deconvolution workflows. This matrix maps samples to their respective
#' cell types.
#'
#' The function is particularly useful for preparing reference data for
#' signature matrix generation using CIBERSORTx. Each cell in the output is
#' encoded as:
#' - 0 for 'Unknown' mappings
#' - 1 if the sample belongs to the cell type in that row
#' - 2 if the sample belongs to a different cell type
#'
#' @param data A numeric matrix of pure cell type profiles with genes as row
#'   names and samples as columns. Column names should contain identifiers that
#'   can be mapped to cell types.
#' @param mapping_rules Either: 1. A named list where names are cell type labels
#'   and values are regular expression patterns used to match sample column
#'   names (e.g., list("CD8+ T cells" = "CD8", "Monocytes" = "Mono")), OR 2. A
#'   character vector with the same length and order as colnames(data), directly
#'   specifying the cell type for each sample.
#' @param verbose Logical. If TRUE, displays additional messages during
#'   processing.
#' @param return_format A string specifying the return format: "matrix" or
#'   "tibble" (default).
#'
#' @return If return_format is "matrix", a numeric matrix with cell type groups
#'   as rows and samples as columns. If return_format is "tibble", a tibble with
#'   a "cell_type" column and columns for each sample.
#'
#' @examples
#' \dontrun{
#' # Example using a named list of regex patterns
#' pure_samples <- readRDS(system.file("extdata", "pure_samples_matrix.rds",
#'                                      package = "proteoDeconv"))
#'
#' mapping_rules <- list(
#'   "CD8+ T cells" = "CD8",
#'   "Monocytes" = "Mono"
#' )
#'
#' phenoclasses1 <- create_phenoclasses(
#'   data = pure_samples,
#'   mapping_rules = mapping_rules,
#'   verbose = TRUE
#' )
#'
#' # Example using a character vector of direct cell type assignments
#' cell_types <- c("CD8+ T cells", "CD8+ T cells", "Monocytes", "Monocytes", "Unknown")
#' phenoclasses2 <- create_phenoclasses(
#'   data = pure_samples,
#'   mapping_rules = cell_types,
#'   verbose = TRUE
#' )
#' }
#'
#' @seealso \code{\link{create_signature_matrix}} which uses the phenoclasses
#'   matrix to generate a cell type signature matrix.
#'
#' @export
create_phenoclasses <- function(
  data,
  mapping_rules,
  verbose = FALSE,
  return_format = c("tibble", "matrix")
) {
  return_format <- match.arg(return_format)

  if (!is.matrix(data)) {
    stop("data must be a matrix")
  }

  cols <- colnames(data)
  if (is.null(cols)) {
    stop("Matrix must have column names representing sample IDs")
  }

  if (is.list(mapping_rules)) {
    if (verbose) {
      message("Using regex-based mapping rules")
    }
    group_mapped <- map_cell_groups(
      cols,
      mapping_rules,
      default_group = "Unknown",
      verbose = verbose
    )
    cell_type_to_group <- stats::setNames(group_mapped, cols)
  } else if (is.character(mapping_rules)) {
    if (length(mapping_rules) != length(cols)) {
      stop(
        "When providing a character vector, its length must match the number of columns in the data matrix"
      )
    }
    if (verbose) {
      message("Using direct cell type assignments")
    }
    cell_type_to_group <- stats::setNames(mapping_rules, cols)
  } else {
    stop(
      "mapping_rules must be either a named list of regex patterns or a character vector of cell type assignments"
    )
  }

  valid_cols <- names(cell_type_to_group)
  valid_groups <- unique(cell_type_to_group)

  valid_groups <- setdiff(valid_groups, "Unknown")

  if (length(valid_groups) == 0) {
    stop("No valid cell types found after mapping. Check your mapping_rules.")
  }

  phenotype_classes <- matrix(
    NA,
    nrow = length(valid_groups),
    ncol = length(valid_cols),
    dimnames = list(valid_groups, valid_cols)
  )

  for (group in valid_groups) {
    for (j in seq_along(valid_cols)) {
      current_group <- cell_type_to_group[[valid_cols[j]]]
      if (current_group == "Unknown") {
        phenotype_classes[group, j] <- 0
      } else if (current_group == group) {
        phenotype_classes[group, j] <- 1
      } else {
        phenotype_classes[group, j] <- 2
      }
    }
  }

  if (return_format == "matrix") {
    return(phenotype_classes)
  } else {
    phenotype_classes_df <- as.data.frame(phenotype_classes) |>
      tibble::rownames_to_column(var = "cell_type") |>
      tibble::as_tibble()
    return(phenotype_classes_df)
  }
}
