#' Signature matrix based on data from Rieckmann et al. (2017)
#'
#' This is one of the signature matrices used in the proteoDeconv manuscript
#' (Zamore et al., 2025), based on proteomics data from Rieckmann et al. (2017).
#'
#' @format An `.rds` file containing a matrix (proteins as rows,
#' immune cell types as columns).
#'
#' @examples
#' path <- system.file("extdata", "signature_rieckmann.rds", package = "proteoDeconv")
#' sig <- readRDS(path)
"signature_rieckmann"
