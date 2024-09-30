#' Gene Symbols Map
#'
#' Data of approved gene symbols.
#'
#' @format ## `gene_symbols_map`
#' A data frame with 103,478 rows and 2 columns:
#' \describe{
#'   \item{Symbol}{Gene symbol}
#'   \item{Approved.Symbol}{Approved gene symbol}
#' }
#' @source HGNC

gene_symbols_map <- HGNChelper::getCurrentHumanMap()

usethis::use_data(gene_symbols_map, overwrite = TRUE)
