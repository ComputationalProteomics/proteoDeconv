#' Get Signature Markers
#'
#' Extracts unique signature markers from a signature matrix.
#'
#' @param Y A list of data matrices.
#' @param signature_matrix A matrix of signature marker values.
#'
#' @return A character matrix of markers.
#' @keywords internal
get_signature_markers <- function(Y, signature_matrix) {
    if (length(Y) > 1) {
        genes <- unique(c(rownames(Y[[1]]), rownames(Y[[2]])))
    } else {
        genes <- rownames(Y[[1]])
    }
    cell_types <- colnames(signature_matrix)
    index_matrix <- NULL
    for (s in seq_len(ncol(signature_matrix))) {
        for (k in seq_len(ncol(signature_matrix))) {
            if (k != s) {
                i <- (signature_matrix[, s] > 1000 & signature_matrix[, s] > 3 * signature_matrix[, k])
                marker_unique <- rownames(signature_matrix)[i]
                marker_unique <- marker_unique[!is.na(match(marker_unique, genes))]
                if (length(marker_unique) >= 1) {
                    new_rows <- cbind(
                        rep(cell_types[s], length(marker_unique)),
                        rep(cell_types[k], length(marker_unique)),
                        marker_unique
                    )
                    index_matrix <- if (is.null(index_matrix)) {
                        new_rows
                    } else {
                        rbind(index_matrix, new_rows)
                    }
                }
            }
        }
    }
    return(index_matrix)
}

#' Transform BayesDeBulk Output
#'
#' Transforms BayesDeBulk cell fraction output into a tidy long data frame.
#'
#' @param bayes_cell_fractions A matrix of cell fractions from BayesDeBulk.
#' @param ... Additional arguments for \code{dplyr::mutate}.
#'
#' @return A tibble in long format with columns for sample, cell type, cell count, and method.
#' @keywords internal
transform_bayesdebulk_output <- function(bayes_cell_fractions, ...) {
    df <- as.data.frame(bayes_cell_fractions) |>
        tibble::rownames_to_column("sample")
    df_long <- tidyr::pivot_longer(
        df,
        cols = -sample,
        names_to = "cell_type",
        values_to = "cell_count"
    ) |>
        dplyr::mutate(method = "BayesDeBulk") |>
        dplyr::mutate(...)
    return(df_long)
}

#' Deconvolute Bulk Proteome Data using BayesDeBulk
#'
#' Deconvolutes bulk proteome data using the BayesDeBulk algorithm.
#'
#' @param data A matrix or data frame of bulk proteome data.
#' @param signature_matrix A matrix containing signature marker values.
#' @param n_iter Number of iterations for the deconvolution; default is 1000.
#' @param burn_in Number of burn-in iterations; default is 100.
#' @param ... Additional arguments passed to \code{dplyr::mutate} when transforming the output.
#'
#' @return A tibble in long format containing cell type fractions, sample identifiers, and additional annotations.
#'
#' @details This function calculates signature markers using \code{get_signature_markers}, runs the
#' BayesDeBulk deconvolution, and then transforms the output using \code{transform_bayesdebulk_output}.
#'
#' @examples
#' \dontrun{
#' result <- deconvolute_bayesdebulk(bulk_data, signature_matrix)
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @export
deconvolute_bayesdebulk <- function(data,
                                    signature_matrix,
                                    n_iter = 1000,
                                    burn_in = 100,
                                    ...) {
    if (!requireNamespace("BayesDeBulk", quietly = TRUE)) {
      stop("Package 'BayesDeBulk' is required but not installed.")
    }
    markers <- get_signature_markers(list(data), signature_matrix)
    bayes <- BayesDeBulk::BayesDeBulk(
        n.iter = n_iter,
        burn.in = burn_in,
        Y = list(data),
        markers = markers
    )
    bayes_cell_fractions <- bayes$cell.fraction
    result <- transform_bayesdebulk_output(bayes_cell_fractions, ...)
    return(result)
}
