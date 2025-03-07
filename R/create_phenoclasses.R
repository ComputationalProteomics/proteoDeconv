#' Create phenoclasses matrix for cell type deconvolution
#'
#' Generates a phenotype classification matrix from sample data for use in cell
#' type deconvolution workflows. This matrix maps samples to their respective
#' cell types using a set of mapping rules.
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
#' @param mapping_rules A named list where names are cell type labels and values
#'   are regular expression patterns used to match sample column names (e.g.,
#'   list("CD8+ T cells" = "CD8", "Monocytes" = "Mono")).
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
#' pure_samples <- readRDS(system.file("extdata", "pure_samples_matrix.rds",
#'                                      package = "proteoDeconv"))
#'
#' mapping_rules <- list(
#'   "CD8+ T cells" = "CD8",
#'   "Monocytes" = "Mono"
#' )
#'
#' phenoclasses <- create_phenoclasses(
#'   data = pure_samples,
#'   mapping_rules = mapping_rules,
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
    stop("immune_cells_mat must be a matrix")
  }

  cols <- colnames(data)
  if (is.null(cols)) {
    stop("Matrix must have column names representing sample IDs")
  }

  group_mapped <- map_cell_groups(
    cols,
    mapping_rules,
    default_group = "Unknown",
    verbose = verbose
  )
  cell_type_to_group <- stats::setNames(group_mapped, cols)

  valid_cols <- names(cell_type_to_group)
  valid_groups <- unique(cell_type_to_group)

  valid_groups <- setdiff(valid_groups, "Unknown")

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
