#' Map Cell Groups from Column Names
#'
#' Maps a character vector of column names to predefined cell group categories using regex patterns.
#'
#' @param column_names A character vector of column names.
#' @param mapping_rules A named list of regex patterns.
#' @param default_group A string that represents the default group name for columns not matching any pattern.
#' @param verbose Logical. If TRUE, prints mapping process messages.
#'
#' @return A character vector with the mapped cell group for each column name.
#'
#' @examples
#' column_names <- c("Bcell_1", "T4.naive", "unknown")
#' map_cell_groups(column_names, verbose = TRUE)
#' @export

map_cell_groups <- function(column_names, mapping_rules, default_group = "Unknown", verbose = FALSE) {
    if (!is.character(column_names)) {
        stop("`column_names` must be a character vector.")
    }
    n <- length(column_names)

    final_mapping <- rep(default_group, n)
    for (group in names(mapping_rules)) {
        patterns <- mapping_rules[[group]]
        combined_pattern <- paste(patterns, collapse = "|")
        matches <- stringr::str_detect(column_names, stringr::regex(combined_pattern, ignore_case = TRUE))
        idx <- which(matches & (final_mapping == default_group))
        final_mapping[idx] <- group

        if (verbose && length(idx) > 0) {
            message(sprintf("Mapping to '%s': %s", group, paste(column_names[idx], collapse = ", ")))
        }
    }

    if (verbose) {
        unmatched_idx <- which(final_mapping == default_group)
        if (length(unmatched_idx) > 0) {
            message("Unmapped columns: ", paste(column_names[unmatched_idx], collapse = ", "))
        }
    }

    return(final_mapping)
}
