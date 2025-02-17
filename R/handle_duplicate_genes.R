#' Handle Duplicate Genes
#'
#' Resolves duplicate gene entries in an expression data frame using the chosen mode.
#'
#' @param data A data frame containing gene expression data.
#' @param gene_column A string specifying the name of the gene identifier column. Default is "Genes".
#' @param duplicate_mode A string specifying the approach to handle duplicates. Options are
#'   \code{"slice"} (keep the row with the maximum median value) or \code{"merge"} (merge rows by taking the median of numeric columns).
#'
#' @return A data frame with duplicate gene rows handled.
#'
#' @details For \code{slice} mode the function calculates the median value across all numeric columns per row,
#' then retains only the row with the maximum median for duplicate gene identifiers.
#' In \code{merge} mode, duplicates are grouped by the gene identifier and numeric columns are summarized by their median.
#'
#' @importFrom dplyr rowwise mutate slice_max distinct filter ungroup group_by summarise across
#' @importFrom stats median
#' @importFrom glue glue
#' @export
handle_duplicate_genes <- function(data, gene_column = "Genes", duplicate_mode) {
    data <- handle_input_data(data, gene_column)
    if (!gene_column %in% colnames(data)) {
        stop(glue::glue("Column '{gene_column}' not found in the data."))
    }

    data <- switch(duplicate_mode,
        slice = {
            data |>
                dplyr::rowwise() |>
                dplyr::mutate(median = stats::median(dplyr::c_across(tidyselect::where(is.numeric)), na.rm = TRUE)) |>
                dplyr::slice_max(median) |>
                dplyr::distinct(.data[[gene_column]], .keep_all = TRUE) |>
                dplyr::filter(!is.na(.data[[gene_column]])) |>
                dplyr::select(-median) |>
                dplyr::ungroup()
        },
        merge = {
            data |>
                dplyr::filter(!is.na(.data[[gene_column]])) |>
                dplyr::group_by(.data[[gene_column]]) |>
                dplyr::summarise(across(where(is.numeric), median, na.rm = TRUE)) |>
                dplyr::ungroup()
        },
        stop(glue::glue("Unsupported duplicate mode: {duplicate_mode}"))
    )

    return(data)
}
