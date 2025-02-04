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
    data <- proteoDeconv::handle_input_data(data, gene_column)

    if (!gene_column %in% colnames(data)) {
        stop(glue::glue("Column '{gene_column}' not found in the data."))
    }

    numeric_data <- data |> dplyr::select(-{{ gene_column }})

    out_data <- switch(imputation_mode,
        "lowest_value" = {
            lowest_value <- min_nonzero(numeric_data)
            numeric_data <- replace(numeric_data, is.na(numeric_data), lowest_value)
            numeric_data |> dplyr::mutate({{ gene_column }} := data[[gene_column]])
        },
        "knn" = impute_helper(numeric_data, "knn", gene_column, data),
        "zero" = impute_helper(numeric_data, "zero", gene_column, data),
        "MLE" = impute_helper(numeric_data, "MLE", gene_column, data),
        "bpca" = impute_helper(numeric_data, "bpca", gene_column, data),
        "RF" = impute_helper(numeric_data, "RF", gene_column, data),
        "min" = impute_helper(numeric_data, "min", gene_column, data),
        "MinDet" = impute_helper(numeric_data, "MinDet", gene_column, data),
        "MinProb" = impute_helper(numeric_data, "MinProb", gene_column, data),
        "QRILC" = impute_helper(numeric_data, "QRILC", gene_column, data),
        "mixed" = impute_helper(numeric_data, "mixed", gene_column, data),
        "nbavg" = impute_helper(numeric_data, "nbavg", gene_column, data),
        "with" = impute_helper(numeric_data, "with", gene_column, data, val = 1),
        stop(glue::glue("Unsupported imputation mode: {imputation_mode}"))
    )


    out_data
}

impute_helper <- function(numeric_data, method, gene_column, data, ...) {
    out <- MsCoreUtils::impute_matrix(as.matrix(numeric_data), method = method, ...)
    tibble::as_tibble(out) |>
        dplyr::mutate({{ gene_column }} := data[[gene_column]], .before = 1)
}
