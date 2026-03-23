# Simulate artifical mixtures from bulk proteome measurements

This function simulates bulk proteome data by mixing bulk sample
measurements.

## Usage

``` r
simulate_data(
  data,
  cell_types,
  seed = NULL,
  ncells = 100,
  nsamples = 100,
  filter_genes = TRUE,
  scenario = "random",
  whitelist = NULL,
  blacklist = NULL
)
```

## Arguments

- data:

  A numeric matrix containing protein abundance data with protein
  identifiers as row names and samples as columns.

- cell_types:

  A character vector indicating the cell type associated with each
  column in the input matrix. Must have the same length as the number of
  columns in `data`.

- seed:

  An integer used as random seed for reproducibility. Default is NULL
  (no seed).

- ncells:

  Integer specifying the number of cells to use for each simulated bulk
  sample.

- nsamples:

  Integer specifying the number of bulk samples to simulate. Default is
  100.

- filter_genes:

  Logical; whether to filter out proteins/genes with low expression
  before simulation. Default is TRUE.

- scenario:

  String specifying the simulation scenario:

  - `"random"` (default): Random cell type proportions for each sample

  - `"even"`: Even proportions of all cell types

  - See SimBu documentation for additional scenarios

- whitelist:

  Optional character vector of cell types to include in the simulation.
  If provided, only these cell types will be used. Default is NULL (use
  all).

- blacklist:

  Optional character vector of cell types to exclude from the
  simulation. Default is NULL (exclude none).

## Value

A list containing two matrices:

- `simulated_data`: Matrix of simulated bulk protein abundance data
  (proteins as rows, simulated samples as columns)

- `cell_fractions`: Matrix of cell type fractions used for each
  simulation (cell types as rows, simulated samples as columns)

## Details

This function uses the [SimBu package](https://omnideconv.org/SimBu/) to
generate synthetic bulk samples by artifically mixing samples.

## Examples

``` r
# Create example data
cell_data <- matrix(abs(rnorm(1500, mean = 500, sd = 200)), nrow = 100, ncol = 15)
rownames(cell_data) <- paste0("Protein", 1:100)
colnames(cell_data) <- paste0("Cell", 1:15)

# Define cell types
cell_types <- rep(c("T_cell", "B_cell", "Monocyte"), each = 5)

# Run simulation
if (FALSE) { # \dontrun{
sim_results <- simulate_data(
  data = cell_data,
  cell_types = cell_types,
  seed = 42,
  nsamples = 20,
  scenario = "random"
)

dim(sim_results$simulated_data)
dim(sim_results$cell_fractions)
} # }
```
