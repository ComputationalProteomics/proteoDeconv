#' Update Gene Symbols
#'
#' This function updates gene symbols in a given dataset using a predefined gene symbols map.
#'
#' @param data A data frame or vector containing gene symbols. If a data frame is provided, the first column is assumed to contain the gene symbols.
#' @param gene_column A string specifying the name of the column containing gene symbols. Default is NULL. If NULL, the first column is used.
#' @param filepath A string specifying the file path for verbose output. Default is an empty string.
#' @param verbose A logical value indicating whether to print detailed messages. Default is FALSE.
#'
#' @return If `data` is a data frame, the function returns a data frame with updated gene symbols. If `data` is a vector, the function returns a vector of updated gene symbols.
#'
#' @details The function uses the `HGNChelper` package to check and update gene symbols based on a predefined gene symbols map. If `verbose` is TRUE, the function prints detailed messages about the number of approved, not approved, and NA gene symbols.
#'
#' @importFrom dplyr rename pull mutate filter
#' @importFrom glue glue
#' @importFrom HGNChelper checkGeneSymbols
#' @export
update_gene_symbols <- function(data, gene_column = NULL, filepath = "", verbose = FALSE) {
  if (!is.vector(data)) {
    if (is.null(gene_column)) {
      gene_column <- colnames(data)[1]
    }
    data <- handle_input_data(data, gene_column)
    checked <- HGNChelper::checkGeneSymbols(data[[gene_column]], map = gene_symbols_map)
  } else {
    checked <- HGNChelper::checkGeneSymbols(data, map = gene_symbols_map)
  }

  if (verbose) {
    message(glue::glue(
      "File: {filepath}\n",
      "Number of approved symbols: {sum(checked$Approved, na.rm = TRUE)}\n",
      "Number of not approved symbols: {sum(!checked$Approved, na.rm = TRUE)}\n",
      "Number of NA in suggested symbols: {sum(is.na(checked$Suggested.Symbol))}"
    ))
  }

  if (!is.vector(data)) {
    data <- data |>
      dplyr::mutate(!!gene_column := checked$Suggested.Symbol) |>
      dplyr::filter(!is.na(!!sym(gene_column)))
    return(data)
  }

  return(checked$Suggested.Symbol)
}
