#' Simulate artifical mixtures from bulk proteome measurements
#'
#' This function simulates bulk proteome data by mixing bulk sample measurements.
#'
#' @param data A numeric matrix containing protein abundance data with protein identifiers
#'        as row names and samples as columns.
#' @param cell_types A character vector indicating the cell type associated with each column
#'        in the input matrix. Must have the same length as the number of columns in \code{data}.
#' @param seed An integer used as random seed for reproducibility. Default is NULL (no seed).
#' @param ncells Integer specifying the number of cells to use for each simulated bulk sample.
#' @param nsamples Integer specifying the number of bulk samples to simulate. Default is 100.
#' @param filter_genes Logical; whether to filter out proteins/genes with low expression
#'        before simulation. Default is TRUE.
#' @param scenario String specifying the simulation scenario:
#'        \itemize{
#'          \item \code{"random"} (default): Random cell type proportions for each sample
#'          \item \code{"even"}: Even proportions of all cell types
#'          \item See SimBu documentation for additional scenarios
#'        }
#' @param whitelist Optional character vector of cell types to include in the simulation.
#'        If provided, only these cell types will be used. Default is NULL (use all).
#' @param blacklist Optional character vector of cell types to exclude from the simulation.
#'        Default is NULL (exclude none).
#'
#' @return A list containing two matrices:
#'   \itemize{
#'     \item \code{simulated_data}: Matrix of simulated bulk protein abundance data
#'           (proteins as rows, simulated samples as columns)
#'     \item \code{cell_fractions}: Matrix of cell type fractions used for each simulation
#'           (cell types as rows, simulated samples as columns)
#'   }
#'
#' @details This function uses the [SimBu package](https://omnideconv.org/SimBu/) to generate synthetic bulk
#' samples by artifically mixing samples.
#'
#'
#' @examples
#' # Create example data
#' cell_data <- matrix(abs(rnorm(1500, mean = 500, sd = 200)), nrow = 100, ncol = 15)
#' rownames(cell_data) <- paste0("Protein", 1:100)
#' colnames(cell_data) <- paste0("Cell", 1:15)
#'
#' # Define cell types
#' cell_types <- rep(c("T_cell", "B_cell", "Monocyte"), each = 5)
#'
#' # Run simulation
#' \dontrun{
#' sim_results <- simulate_data(
#'   data = cell_data,
#'   cell_types = cell_types,
#'   seed = 42,
#'   nsamples = 20,
#'   scenario = "random"
#' )
#'
#' dim(sim_results$simulated_data)
#' dim(sim_results$cell_fractions)
#' }
#' @export
simulate_data <- function(
  data,
  cell_types,
  seed = NULL,
  ncells = 100,
  nsamples = 100,
  filter_genes = TRUE,
  scenario = "random",
  whitelist = NULL,
  blacklist = NULL
) {
  if (!is.matrix(data)) {
    stop("Input must be a matrix")
  }

  if (is.null(rownames(data))) {
    stop("Matrix must have row names as protein identifiers")
  }

  sample_cols <- colnames(data)

  if (length(cell_types) != length(sample_cols)) {
    stop(
      "Length of cell_types must match the number of samples (columns) in the matrix."
    )
  }

  annotation <- data.frame(
    ID = sample_cols,
    cell_type = cell_types,
    stringsAsFactors = FALSE
  )

  counts <- data

  tpm <- handle_scaling(data, unlog = FALSE, tpm = TRUE)

  counts_sparse <- Matrix::Matrix(counts, sparse = TRUE)
  tpm_sparse <- Matrix::Matrix(tpm, sparse = TRUE)

  dataset <- SimBu::dataset(
    annotation = annotation,
    count_matrix = counts_sparse,
    tpm_matrix = tpm_sparse,
    name = "dataset",
    filter_genes = filter_genes
  )

  set.seed(seed)

  simulation <- SimBu::simulate_bulk(
    dataset,
    scenario = scenario,
    scaling_factor = "NONE",
    seed = seed,
    nsamples = nsamples,
    ncells = ncells,
    whitelist = whitelist,
    blacklist = blacklist
  )

  assays_list <- SummarizedExperiment::assays(simulation$bulk)

  simulated_data <- assays_list[["bulk_tpm"]]

  cell_fractions <- simulation$cell_fractions

  list(
    simulated_data = as.matrix(simulated_data),
    cell_fractions = as.matrix(cell_fractions)
  )
}
