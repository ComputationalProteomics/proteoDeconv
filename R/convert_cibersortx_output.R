#' Convert CIBERSORTx Output
#'
#' This function converts the output from CIBERSORTx into a long format.
#'
#' @param data A data frame containing CIBERSORTx output.
#' @param file A character string specifying the file name (optional).
#' @param filter A logical value indicating whether to filter rows with P-value > 0.05.
#' @return A data frame in long format.
#' @importFrom dplyr filter select all_of
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect any_of
#' @importFrom glue glue
#' @export
convert_cibersortx_output <- function(data, file = NULL, filter = FALSE) {
    extra_cols <- c("P.value", "Correlation", "RMSE", "P-value")

    if (!("Mixture" %in% colnames(data))) {
        stop("No Mixture column found")
    }

    if (filter) {
        num_rows <- nrow(data)
        data <- dplyr::filter(data, .data$`P-value` < 0.05)
        diff <- num_rows - nrow(data)
        if (!is.null(file)) {
            message(glue::glue("Removed {diff} samples with CIBERSORTx P-value > 0.05 from {file}"))
        } else {
            message(glue::glue("Removed {diff} samples with CIBERSORTx P-value > 0.05"))
        }
    }

    data |>
        dplyr::select(-tidyselect::any_of(extra_cols)) |>
        tidyr::pivot_longer(
            cols = -dplyr::all_of("Mixture"),
            names_to = "cell_type",
            values_to = "value"
        ) |>
        tidyr::pivot_wider(
            names_from = "Mixture",
            values_from = "value"
        )
}
