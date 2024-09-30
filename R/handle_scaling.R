#' Handle Scaling of Protein Data
#'
#' This function handles the scaling of protein data by optionally unlogging and converting to TPM (Transcripts Per Million).
#'
#' @param data A data frame containing protein data. The first column should contain gene identifiers.
#' @param gene_column A string specifying the name of the column containing gene identifiers. Default is "Genes".
#' @param unlog A logical value indicating whether to unlog the data (convert from log2 scale). Default is TRUE.
#' @param tpm A logical value indicating whether to convert the data to TPM (Transcripts Per Million). Default is TRUE.
#'
#' @return A data frame with the scaled protein data.
#'
#' @importFrom dplyr select mutate across relocate
#' @importFrom glue glue
#' @export
handle_scaling <- function(data, gene_column = "Genes", unlog = TRUE, tpm = TRUE) {
    if (!gene_column %in% colnames(data)) {
        stop(glue::glue("Column '{gene_column}' not found in the data."))
    }

    numeric_data <- data |> dplyr::select(-{{ gene_column }})

    if (unlog) {
        numeric_data <- numeric_data |> dplyr::mutate(across(everything(), ~ 2^.x))
    }

    if (tpm) {
        numeric_data <- numeric_data |> dplyr::mutate(across(where(is.numeric), ~ (.x / sum(.x)) * 1e+06))
    }

    result_data <- numeric_data |>
        dplyr::mutate({{ gene_column }} := data[[gene_column]]) |>
        dplyr::relocate({{ gene_column }}, .before = 1)

    return(result_data)
}
