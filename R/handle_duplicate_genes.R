#' Handle Duplicate Genes
#'
#' This function handles duplicate genes in a data frame using specified modes.
#'
#' @param data A data frame containing gene expression data.
#' @param gene_column A string specifying the name of the column containing gene identifiers. Default is "Genes".
#' @param duplicate_mode A string specifying the mode to handle duplicates. Options are "slice" and "merge".
#'
#' @return A data frame with duplicates handled according to the specified mode.
#'
#' @details The function handles duplicate genes by either slicing to keep the row with the maximum median value or merging duplicates by taking the median of numeric columns.
#'
#' @importFrom dplyr rowwise mutate slice_max distinct filter ungroup group_by summarise_all across
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
            data <- data |>
                dplyr::filter(!is.na(.data[[gene_column]])) |>
                dplyr::group_by(.data[[gene_column]]) |>
                dplyr::summarise(across(where(is.numeric), median, na.rm = TRUE)) |>
                dplyr::ungroup()
        },
        stop(glue::glue("Unsupported duplicate mode: {duplicate_mode}"))
    )

    return(data)
}
