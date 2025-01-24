#' Handle Input Data
#'
#' This function handles various input data types and normalizes them into a common format with a "Genes" column.
#'
#' @param x An input data object. Supported types are SummarizedExperiment, tibble, or data frame.
#' @param gene_column A string specifying the name of the column containing gene identifiers. Default is "Genes".
#' @param as_tibble A logical value indicating whether to return a tibble with the first column being Genes (TRUE) or a data frame with the Genes column as rownames (FALSE). Default is TRUE.
#'
#' @return A data frame or tibble with the input data normalized to have a "Genes" column.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column column_to_rownames
#' @export
handle_input_data <- function(x, gene_column = "Genes", as_tibble = TRUE) {
    if (inherits(x, "SummarizedExperiment")) {
        df <- as.data.frame(SummarizedExperiment::assay(x))
        if (!is.null(rownames(df))) {
            df <- tibble::rownames_to_column(df, var = gene_column)
        }
        if (!as_tibble) {
            df <- tibble::column_to_rownames(df, var = gene_column)
        }
        return(df)
    }
    if (is.data.frame(x)) {
        # If rownames are present, turn them into a column
        if (!gene_column %in% colnames(x) && !is.null(rownames(x))) {
            x <- tibble::rownames_to_column(x, var = gene_column)
        }
        if (!as_tibble) {
            x <- tibble::column_to_rownames(x, var = gene_column)
        }
        return(x)
    }
    stop("Unsupported input type.")
}
