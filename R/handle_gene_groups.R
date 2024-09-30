#' Handle Protein/Gene Groups
#'
#' This function handles protein/gene groups by picking the first protein/gene in each group.
#'
#' @param data A data frame containing protein data. The specified gene column should contain protein/gene groups separated by a specified separator.
#' @param gene_column A string specifying the name of the column containing gene groups. Default is "Genes".
#' @param separator A string specifying the separator used to split protein/gene groups. Default is ";".
#'
#' @return A data frame with the first protein/gene in each group.
#'
#' @details The function splits the protein/gene groups in the specified column using the provided separator and picks the first protein/gene in each group.
#'
#' @importFrom glue glue
#' @export
handle_gene_groups <- function(data, gene_column = "Genes", separator = ";") {
    if (!gene_column %in% colnames(data)) {
        stop(glue::glue("Column '{gene_column}' not found in the data."))
    }
    data[[gene_column]] <- strsplit(data[[gene_column]], separator) |>
        lapply(function(x) x[1]) |>
        unlist()
    data
}
