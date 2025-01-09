#' Handle Missing Values in a Data Frame
#'
#' This function handles missing values in a data frame by either replacing them with the lowest non-zero value
#' or using k-nearest neighbors (KNN) imputation.
#'
#' @param data A data frame containing the data to be processed.
#' @param gene_column A string specifying the name of the column containing gene identifiers. Default is "Genes".
#' @param imputation_mode A string specifying the imputation method to use. Options are "lowest_value" or "knn". Default is "lowest_value".
#'
#' @return A data frame with missing values imputed.
#'
#' @importFrom dplyr select mutate relocate
#' @importFrom glue glue
#' @importFrom MsCoreUtils impute_matrix
#' @importFrom tibble as_tibble
#' @export
handle_missing_values <- function(data, gene_column = "Genes", imputation_mode = "lowest_value") {
    if (!gene_column %in% colnames(data)) {
        stop(glue::glue("Column '{gene_column}' not found in the data."))
    }

    numeric_data <- data |> dplyr::select(-{{ gene_column }})

    data <- switch(imputation_mode,
        lowest_value = {
            lowest_value <- min_nonzero(numeric_data)
            numeric_data <- replace(numeric_data, is.na(numeric_data), lowest_value)
            numeric_data |> dplyr::mutate({{ gene_column }} := data[[gene_column]])
        },
        median = {
            median_value <- median(as.matrix(numeric_data), na.rm = TRUE)
            numeric_data <- replace(numeric_data, is.na(numeric_data), median_value)
            numeric_data |> dplyr::mutate({{ gene_column }} := data[[gene_column]])
        },
        knn = {
            m <- as.matrix(numeric_data)
            m <- MsCoreUtils::impute_matrix(m, method = "knn") |> tibble::as_tibble()
            m |> dplyr::mutate({{ gene_column }} := data[[gene_column]])
        },
        zero = {
            numeric_data[is.na(numeric_data)] <- 0
            numeric_data |> dplyr::mutate({{ gene_column }} := data[[gene_column]])
        },
        stop(glue::glue("Unsupported imputation mode: {imputation_mode}"))
    )

    data |> dplyr::relocate({{ gene_column }}, .before = 1)
}

min_nonzero <- function(x) {
    non_zero_values <- x[x != 0]
    min(non_zero_values, na.rm = TRUE)
}
