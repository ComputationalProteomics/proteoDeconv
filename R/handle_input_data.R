#' Handle input data
#'
#' This function processes input data, automatically identifying the gene column
#' if not specified. It supports matrices, data frames, and SummarizedExperiment objects.
#'
#' @param x A matrix, data frame, or SummarizedExperiment object.
#' @param gene_column Optional. Name of the column containing gene identifiers.
#' @param as_tibble Logical. If FALSE, moves the gene column to row names.
#'
#' @return A processed data frame or tibble.
#' @export
handle_input_data <- function(x, gene_column = NULL, as_tibble = TRUE) {
    guess_gene_col <- function(df) {
        char_cols <- sapply(df, is.character)
        if (!any(char_cols)) {
            return(NULL)
        }
        names(df)[which(char_cols)[1]]
    }

    if (inherits(x, "SummarizedExperiment")) {
        df <- as.data.frame(SummarizedExperiment::assay(x))
    } else if (is.matrix(x)) {
        df <- as.data.frame(x, stringsAsFactors = FALSE)
    } else if (is.data.frame(x)) {
        df <- x
    } else {
        stop("Unsupported input type. Provide a matrix, data.frame, or SummarizedExperiment.")
    }

    if (!is.null(gene_column)) {
        if (!gene_column %in% names(df)) {
            guessed_col <- guess_gene_col(df)
            if (!is.null(guessed_col)) {
                colnames(df)[colnames(df) == guessed_col] <- gene_column
            } else if (!is.null(rownames(df)) && any(nzchar(rownames(df)))) {
                df <- tibble::rownames_to_column(df, var = gene_column)
            } else {
                message("No suitable gene column found.")
            }
        }
    } else {
        guessed_col <- guess_gene_col(df)
        if (!is.null(guessed_col)) {
            colnames(df)[colnames(df) == guessed_col] <- "Genes"
            gene_column <- "Genes"
        } else if (!is.null(rownames(df)) && any(nzchar(rownames(df)))) {
            df <- tibble::rownames_to_column(df, var = "Genes")
            gene_column <- "Genes"
        } else {
            message("Could not identify a gene column.")
        }
    }

    if (!as_tibble && !is.null(gene_column) && gene_column %in% names(df)) {
        df <- tibble::column_to_rownames(df, var = gene_column)
    }

    df
}
