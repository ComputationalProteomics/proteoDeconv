#' Deconvolute Bulk Data using EPIC
#'
#' Runs the EPIC algorithm on preprocessed proteomic data to estimate cell type fractions.
#'
#' @param preprocessed_data A matrix or data frame of the bulk proteome data.
#' @param signature_df A matrix or data frame of reference signature profiles.
#' @param with_other_cells Logical; if TRUE EPIC will include other cell types. Default is TRUE.
#' @param method_label A character string label used to tag the output with the deconvolution method.
#'
#' @return A tibble in long format with columns for \code{sample}, \code{cell_type}, \code{cell_count}, and \code{method} (with additional fields as provided in \code{...}).
#'
#' @details The function preprocesses both the bulk data and the signature matrix using \code{handle_scaling}
#' and \code{handle_input_data}. It then runs the EPIC deconvolution method from the EPIC package and reshapes
#' the resulting mRNA proportions into a tidy format.
#'
#' @examples
#' \dontrun{
#' result <- deconvolute_epic(bulk_data, signature_matrix, with_other_cells = TRUE, method_label = "EPIC")
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @export
deconvolute_epic <- function(preprocessed_data, signature_df,
                             with_other_cells = TRUE, method_label) {
    signature_df_epic <- handle_scaling(signature_df, tpm = TRUE, unlog = FALSE)
    preprocessed_data_epic <- handle_scaling(preprocessed_data, tpm = TRUE, unlog = FALSE)

    preprocessed_data <- proteoDeconv::handle_input_data(preprocessed_data, as_tibble = FALSE)
    signature_df <- proteoDeconv::handle_input_data(signature_df, as_tibble = FALSE)
    preprocessed_data_epic <- proteoDeconv::handle_input_data(preprocessed_data_epic, as_tibble = FALSE)
    signature_df_epic <- proteoDeconv::handle_input_data(signature_df_epic, as_tibble = FALSE)

    epic_signature <- list(
        refProfiles = as.matrix(signature_df_epic),
        sigGenes = rownames(signature_df_epic)
    )

    epic_res <- EPIC::EPIC(
        preprocessed_data_epic,
        withOtherCells = with_other_cells,
        reference = epic_signature
    )
    result <- epic_res$mRNAProportions %>%
        tibble::as_tibble(rownames = "sample") %>%
        tidyr::pivot_longer(
            cols = -sample,
            names_to = "cell_type",
            values_to = "cell_count"
        ) %>%
        dplyr::mutate(method = method_label)

    return(result)
}
