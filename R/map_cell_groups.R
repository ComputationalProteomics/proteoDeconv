#' Map column names to cell type groups using patterns
#'
#' Maps a character vector of column names to predefined cell type categories by
#' matching against patterns, with optional regex support.
#'
#' @param column_names A character vector of column names to be mapped to cell
#'   type groups.
#' @param mapping_rules A named list where names are the target cell type groups
#'   and values are character vectors of patterns that identify columns
#'   belonging to each group.
#' @param default_group A string that represents the default group name for
#'   columns not matching any pattern. Default is "Unknown".
#' @param verbose Logical. If TRUE, prints mapping process messages showing
#'   which columns were mapped to which groups and which columns remained
#'   unmapped. Default is FALSE.
#' @param use_regex Logical. If TRUE (default), treats patterns in mapping_rules
#'   as regular expressions. If FALSE, performs exact string matching.
#'
#' @return A character vector with the same length as \code{column_names},
#'   containing the mapped cell type group for each column name. Columns that
#'   don't match any pattern will be assigned the \code{default_group}.
#'
#' @details This function iterates through the \code{mapping_rules} list and
#'   attempts to match each column name against the patterns for each cell type
#'   group. The first matching group in the order of the list will be assigned.
#'
#'   When \code{use_regex = TRUE}, patterns are treated as regular expressions
#'   and matching is case-insensitive. When \code{use_regex = FALSE}, exact
#'   string matching is performed, which is case-sensitive.
#'
#'
#' @examples
#' # Example column names from a dataset
#' cols <- c("CD8_T_cell", "T_cell", "B_cell",
#'           "NK_cell", "Monocyte", "Unknown_cell")
#'
#' # Define simple mapping rules for cell types
#' mapping <- list(
#'   "T_cell" = c("CD8_T_cell", "T_cell"),
#'   "B_cell" = c("B_cell"),
#'   "NK_cell" = c("NK_cell"),
#'   "Monocyte" = c("Monocyte")
#' )
#'
#' # Map column names to cell types using exact matching
#' cell_types <- map_cell_groups(cols, mapping, use_regex = FALSE)
#' print(data.frame(column = cols, cell_type = cell_types))
#'
#' # Define mapping rules using regex patterns
#' mapping_regex <- list(
#'   "T_cell" = c("CD8.*", "T_cell"),
#'   "B_cell" = c("B[-_]?cell"),
#'   "NK_cell" = c("NK[-_]?cell"),
#'   "Monocyte" = c("Mono.*")
#' )
#'
#' # Map using regex patterns (default)
#' cell_types_regex <- map_cell_groups(cols, mapping_regex)
#' print(data.frame(column = cols, cell_type = cell_types_regex))
#'
#' @export
map_cell_groups <- function(
  column_names,
  mapping_rules,
  default_group = "Unknown",
  verbose = FALSE,
  use_regex = TRUE
) {
  if (!is.character(column_names)) {
    stop("`column_names` must be a character vector.")
  }
  n <- length(column_names)

  final_mapping <- rep(default_group, n)

  for (group in names(mapping_rules)) {
    patterns <- mapping_rules[[group]]

    if (use_regex) {
      combined_pattern <- paste(patterns, collapse = "|")
      matches <- stringr::str_detect(
        column_names,
        stringr::regex(combined_pattern, ignore_case = TRUE)
      )
    } else {
      matches <- sapply(column_names, function(col) {
        any(col %in% patterns)
      })
    }

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
