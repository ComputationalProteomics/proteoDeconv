#' Update gene symbols in a matrix to approved HGNC nomenclature
#'
#' This function checks and updates gene symbols in a matrix's row names to ensure they
#' conform to the current approved HGNC (HUGO Gene Nomenclature Committee) standards.
#' Non-standard or outdated gene symbols are replaced with their current official symbols.
#'
#' @param data A numeric matrix containing gene expression or protein abundance data
#'        with gene identifiers as row names.
#' @param verbose A logical value indicating whether to print detailed messages about
#'        the number of approved, non-approved, and unmappable gene symbols. Default is FALSE.
#'
#' @return A matrix with updated gene symbols as row names. Rows with unmappable symbols
#'         (where no suggested symbol could be found) are removed from the output.
#'
#' @details The function uses the `HGNChelper` package to check and standardize gene symbols
#' based on a predefined gene symbols map. This is important for
#' ensuring consistency in gene naming across datasets and avoiding issues with outdated
#' or non-standard gene symbols.
#'
#' @export
update_gene_symbols <- function(
  data,
  verbose = FALSE
) {
  if (!is.matrix(data)) {
    stop("Input must be a matrix")
  }

  if (is.null(rownames(data))) {
    stop("Matrix must have row names as gene identifiers")
  }

  gene_symbols <- rownames(data)

  checked <- suppressWarnings(HGNChelper::checkGeneSymbols(
    gene_symbols,
    map = gene_symbols_map
  ))

  if (verbose) {
    message(paste0(
      "Number of approved symbols: ",
      sum(checked$Approved, na.rm = TRUE),
      "\n",
      "Number of not approved symbols: ",
      sum(!checked$Approved, na.rm = TRUE),
      "\n",
      "Number of NA in suggested symbols: ",
      sum(is.na(checked$Suggested.Symbol))
    ))
  }

  valid_symbols <- !is.na(checked$Suggested.Symbol)
  filtered_mat <- data[valid_symbols, , drop = FALSE]

  rownames(filtered_mat) <- checked$Suggested.Symbol[valid_symbols]

  return(filtered_mat)
}
