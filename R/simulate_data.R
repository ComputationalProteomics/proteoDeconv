#' Simulate Bulk Proteome Data
#'
#' This function simulates bulk proteome data based on provided bulk measurements.
#'
#' @param data A data frame containing a "Genes" column and sample columns.
#' @param cell_types A vector indicating the cell type associated with each sample.
#' @param seed Optional seed for reproducibility.
#' @param ncells Number of cells to simulate per sample.
#' @param filter_genes Logical; whether or not to filter genes.
#' @param scenario Simulation scenario, default is "random".
#' @param whitelist Optional whitelist of cell types.
#' @param blacklist Optional blacklist of cell types.
#'
#' @return A list containing simulated proteome data and associated cell fractions.
#' @export
simulate_data <- function(data, cell_types, seed = NULL, ncells = 100,
                          filter_genes = TRUE, scenario = "random",
                          whitelist = NULL, blacklist = NULL) {
    if (!"Genes" %in% names(data)) {
        stop('The input data must contain a "Genes" column.')
    }
    sample_cols <- setdiff(names(data), "Genes")
    if (length(cell_types) != length(sample_cols)) {
        stop("Length of cell_types must match the number of samples (columns excluding 'Genes').")
    }

    annotation <- data.frame(
        ID = sample_cols,
        cell_type = cell_types,
        stringsAsFactors = FALSE
    )

    counts <- data |>
        column_to_rownames("Genes") |>
        as.matrix()

    tpm <- data |>
        handle_scaling(unlog = FALSE, tpm = TRUE) |>
        column_to_rownames("Genes") |>
        as.matrix()

    counts_sparse <- Matrix::Matrix(counts, sparse = TRUE)
    tpm_sparse <- Matrix::Matrix(tpm, sparse = TRUE)

    dataset <- SimBu::dataset(
        annotation = annotation,
        count_matrix = counts_sparse,
        tpm_matrix = tpm_sparse,
        name = "dataset",
        filter_genes = filter_genes
    )

    simulation <- SimBu::simulate_bulk(
        dataset,
        scenario = scenario,
        scaling_factor = "NONE",
        seed = seed,
        nsamples = 100,
        ncells = ncells,
        whitelist = whitelist,
        blacklist = blacklist
    )

    assays_list <- SummarizedExperiment::assays(simulation$bulk)
    if (!"bulk_tpm" %in% names(assays_list)) {
        stop('The assay "bulk_tpm" was not found in the simulation results.')
    }

    simulated_data <- assays_list[["bulk_tpm"]] |>
        as.matrix() |>
        tibble::as_tibble(rownames = "Genes")

    cell_fractions <- simulation$cell_fractions |>
        tibble::as_tibble(rownames = "sample")

    list(
        simulated_data = simulated_data,
        cell_fractions = cell_fractions
    )
}
