#' Map column names to cell type groups using regex patterns
#'
#' Maps a character vector of column names to predefined cell type categories by
#' matching against regex patterns.
#'
#' @param column_names A character vector of column names to be mapped to cell
#'   type groups.
#' @param mapping_rules A named list where names are the target cell type groups
#'   and values are character vectors of regex patterns that identify columns
#'   belonging to each group.
#' @param default_group A string that represents the default group name for
#'   columns not matching any pattern. Default is "Unknown".
#' @param verbose Logical. If TRUE, prints mapping process messages showing
#'   which columns were mapped to which groups and which columns remained
#'   unmapped. Default is FALSE.
#'
#' @return A character vector with the same length as \code{column_names},
#'   containing the mapped cell type group for each column name. Columns that
#'   don't match any pattern will be assigned the \code{default_group}.
#'
#' @details This function iterates through the \code{mapping_rules} list and
#'   attempts to match each column name against the regex patterns for each cell
#'   type group. The first matching group in the order of the list will be
#'   assigned. The function is case-insensitive by default.
#'
#'
#' @examples
#' # Example column names from a dataset
#' cols <- c("CD8_T_cell_donor1", "cd8_tcell_donor2", "NK_cell_sample3",
#'           "b_cell_healthy", "Bcell_patient", "other_cell")
#'
#' # Define mapping rules for cell types
#' mapping <- list(
#'   "T_cell" = c("cd8.*t.*cell", "t.*cell"),
#'   "NK_cell" = c("NK_", "natural.*killer"),
#'   "B_cell" = c("b.*cell")
#' )
#'
#' # Map column names to cell types
#' cell_types <- map_cell_groups(cols, mapping)
#' print(data.frame(column = cols, cell_type = cell_types))
#'
#' @export
map_cell_groups <- function(
  column_names,
  mapping_rules,
  default_group = "Unknown",
  verbose = FALSE
) {
  if (!is.character(column_names)) {
    stop("`column_names` must be a character vector.")
  }
  n <- length(column_names)

  final_mapping <- rep(default_group, n)
  for (group in names(mapping_rules)) {
    patterns <- mapping_rules[[group]]
    combined_pattern <- paste(patterns, collapse = "|")
    matches <- stringr::str_detect(
      column_names,
      stringr::regex(combined_pattern, ignore_case = TRUE)
    )
    idx <- which(matches & (final_mapping == default_group))
    final_mapping[idx] <- group

    if (verbose && length(idx) > 0) {
      message(sprintf(
        "Mapping to '%s': %s",
        group,
        paste(column_names[idx], collapse = ", ")
      ))
    }
  }

  if (verbose) {
    unmatched_idx <- which(final_mapping == default_group)
    if (length(unmatched_idx) > 0) {
      message(
        "Unmapped columns: ",
        paste(column_names[unmatched_idx], collapse = ", ")
      )
    }
  }

  return(final_mapping)
}
