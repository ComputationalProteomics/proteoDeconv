#' Deconvolute bulk proteome data using BayesDeBulk
#'
#' Deconvolutes bulk proteome data using the [BayesDeBulk
#' algorithm](https://github.com/WangLab-MSSM/BayesDeBulk) to estimate cell type
#' proportions in mixed samples based on a signature matrix.
#'
#' @param data A numeric matrix of bulk proteome data with gene identifiers as
#'   row names and samples as columns.
#' @param signature A numeric matrix containing signature marker values with
#'   gene identifiers as row names and cell types as columns.
#' @param n_iter Number of iterations for the MCMC sampling; default is 1000.
#' @param burn_in Number of burn-in iterations to discard; default is 100.
#' @param marker_selection The method to use for marker selection: "limma"
#'   (default), "simple", or a pre-computed marker matrix.
#' @param ... Additional arguments passed to marker selection functions.
#'
#' @return A matrix containing cell type proportions with samples as rows and
#'   cell types as columns.
#'
#' @details This function calculates signature markers using the specified
#'   marker selection method and runs the [BayesDeBulk deconvolution
#'   algorithm](https://github.com/WangLab-MSSM/BayesDeBulk). The marker
#'   selection process identifies genes that are uniquely expressed in specific
#'   cell types, which are then used for the deconvolution.
#'
#' @examples
#' \dontrun{
#' # Load example data and signature matrix
#' data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
#' mixed_samples <- readRDS(data_file)
#'
#' signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
#' signature_matrix <- readRDS(signature_file)
#'
#' # Run deconvolution
#' result <- deconvolute_bayesdebulk(mixed_samples, signature_matrix)
#' }
#'
#' @export
deconvolute_bayesdebulk <- function(
  data,
  signature,
  n_iter = 1000,
  burn_in = 100,
  marker_selection = "limma",
  ...
) {
  if (!requireNamespace("BayesDeBulk", quietly = TRUE)) {
    stop("Package 'BayesDeBulk' is required but not installed.")
  }

  if (!is.matrix(data)) {
    stop("Input 'data' must be a matrix")
  }
  if (!is.matrix(signature)) {
    stop("Input 'signature' must be a matrix")
  }

  if (is.null(rownames(data)) || is.null(rownames(signature))) {
    stop("Both matrices must have gene identifiers as row names")
  }

  data_list <- list(data)

  markers <- get_signature_markers(
    data,
    signature,
    method = marker_selection,
    ...
  )

  bayes <- BayesDeBulk::BayesDeBulk(
    n.iter = n_iter,
    burn.in = burn_in,
    Y = data_list,
    markers = markers
  )

  bayes_cell_fractions <- bayes$cell.fraction

  return(bayes_cell_fractions)
}

#' Get signature markers for deconvolution
#'
#' Extracts unique signature markers from a signature matrix for use in cell type deconvolution.
#' This function selects a method to identify cell type-specific markers.
#'
#' @param Y A numeric matrix of bulk data with gene identifiers as row names.
#' @param signature_mat A numeric matrix of signature marker values with gene identifiers as row names
#'        and cell types as columns.
#' @param method The method to use for marker selection: "limma" (default), "simple",
#'        or a pre-computed marker matrix.
#' @param ... Additional arguments passed to the selected marker selection function.
#'
#' @return A character matrix of markers with columns for cell type, comparison cell type, and gene name.
#' @keywords internal
get_signature_markers <- function(
  Y,
  signature_mat,
  method = "limma",
  ...
) {
  if (is.matrix(method) || is.data.frame(method)) {
    return(as.matrix(method))
  }
  switch(
    method,
    "limma" = get_signature_markers_limma(Y, signature_mat, ...),
    "simple" = get_signature_markers_simple(Y, signature_mat, ...),
    stop(
      "Invalid marker_selection method. Must be 'limma', 'simple', or a pre-computed marker matrix."
    )
  )
}

#' Get signature markers using simple threshold method
#'
#' Extracts unique signature markers from a signature matrix using a simple threshold-based approach.
#' A gene is considered a marker if its expression in one cell type exceeds its expression
#' in other cell types by a specified fold change and meets a minimum expression threshold.
#'
#' @param Y A numeric matrix of bulk data with gene identifiers as row names.
#' @param signature_mat A numeric matrix of signature marker values with gene identifiers as row names
#'        and cell types as columns.
#' @param min_expression Minimum expression level required in the target cell type (default: 1000).
#' @param min_fold_change Minimum fold change required compared to other cell types (default: 3).
#'
#' @return A character matrix of markers with columns for cell type, comparison cell type, and gene name.
#' @keywords internal
get_signature_markers_simple <- function(
  Y,
  signature_mat,
  min_expression = 1000,
  min_fold_change = 3
) {
  if (is.list(Y) && !is.data.frame(Y)) {
    if (length(Y) > 1) {
      genes <- unique(c(rownames(Y[[1]]), rownames(Y[[2]])))
    } else {
      genes <- rownames(Y[[1]])
    }
  } else {
    genes <- rownames(Y)
  }

  cell_types <- colnames(signature_mat)
  index <- NULL

  for (s in seq_len(ncol(signature_mat))) {
    for (k in seq_len(ncol(signature_mat))) {
      if (k != s) {
        i <- (signature_mat[, s] > min_expression &
          signature_mat[, s] > min_fold_change * signature_mat[, k])
        marker_unique <- rownames(signature_mat)[i]
        marker_unique <- marker_unique[!is.na(match(marker_unique, genes))]

        if (length(marker_unique) >= 1) {
          new_rows <- cbind(
            rep(cell_types[s], length(marker_unique)),
            rep(cell_types[k], length(marker_unique)),
            marker_unique
          )
          index <- if (is.null(index)) {
            new_rows
          } else {
            rbind(index, new_rows)
          }
        }
      }
    }
  }

  if (!is.null(index)) {
    colnames(index) <- c("cell_type", "compared_to", "gene")
  }

  return(index)
}

#' Get signature markers using limma-based fold change method
#'
#' Extracts unique signature markers from a signature matrix using a fold-change approach.
#' This method calculates fold changes between cell types and selects genes that have
#' significant differential expression.
#'
#' @param Y A numeric matrix of bulk data with gene identifiers as row names.
#' @param signature_mat A numeric matrix of signature marker values with gene identifiers as row names
#'        and cell types as columns.
#' @param min_fold_change Minimum fold change required for a gene to be considered a marker (default: 3).
#' @param min_expression Minimum expression level required in the target cell type (default: 0).
#' @param pseudo_count Small value added to expression values to avoid division by zero (default: 1e-10).
#'
#' @return A character matrix of markers with columns for cell type, comparison cell type, and gene name.
#' @keywords internal
get_signature_markers_limma <- function(
  Y,
  signature_mat,
  min_fold_change = 3,
  min_expression = 0,
  pseudo_count = 1e-10
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' needed for this function to work. Please install it.")
  }

  genes <- rownames(Y)

  cell_types <- colnames(signature_mat)
  index <- NULL

  for (s in seq_len(ncol(signature_mat))) {
    for (k in seq_len(ncol(signature_mat))) {
      if (k != s) {
        mini <- signature_mat[, c(s, k), drop = FALSE]

        log_fc <- log2((mini[, 1] + 1) / (mini[, 2] + 1))

        fold_change <- mini[, 1] / (mini[, 2] + pseudo_count)

        marker_candidates <- rownames(mini)[
          fold_change > min_fold_change & mini[, 1] > min_expression
        ]

        marker_candidates <- marker_candidates[marker_candidates %in% genes]

        if (length(marker_candidates) > 0) {
          new_rows <- cbind(
            rep(cell_types[s], length(marker_candidates)),
            rep(cell_types[k], length(marker_candidates)),
            marker_candidates
          )

          index <- if (is.null(index)) {
            new_rows
          } else {
            rbind(index, new_rows)
          }
        }
      }
    }
  }

  if (!is.null(index)) {
    colnames(index) <- c("cell_type", "compared_to", "gene")
  }

  return(index)
}
