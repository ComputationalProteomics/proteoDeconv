# Deconvolute protein data using multiple algorithms

This function provides a unified interface for deconvoluting protein or
gene expression data using various algorithms, making it easy to switch
between different deconvolution methods.

## Usage

``` r
deconvolute(algorithm, data, signature = NULL, ...)
```

## Arguments

- algorithm:

  A string specifying the deconvolution algorithm to use. Options are
  "cibersortx", "cibersort", "epic", and "bayesdebulk".

- data:

  A numeric matrix containing protein or gene expression data with
  genes/proteins as row names and samples as columns.

- signature:

  A numeric matrix containing the reference signature with
  genes/proteins as row names and cell types as columns. Required by all
  supported algorithms.

- ...:

  Additional arguments passed to the specific deconvolution functions.

## Value

A numeric matrix with samples as rows and cell types as columns,
representing the estimated proportion of each cell type in each sample.

## Details

The function supports multiple deconvolution algorithms:

- **cibersortx**: Uses Docker-based CIBERSORTx, which requires Docker
  installation and credentials (see
  [`deconvolute_cibersortx`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_cibersortx.md)).

- **cibersort**: Uses the R implementation of CIBERSORT, which requires
  sourcing the original CIBERSORT.R script (see
  [`deconvolute_cibersort`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_cibersort.md)).

- **epic**: Uses the EPIC algorithm which performs constrained least
  squares regression (see
  [`deconvolute_epic`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_epic.md)).

- **bayesdebulk**: Uses the BayesDeBulk Bayesian deconvolution algorithm
  (see
  [`deconvolute_bayesdebulk`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_bayesdebulk.md)).

## See also

[`deconvolute_cibersortx`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_cibersortx.md)
for the CIBERSORTx Docker implementation
[`deconvolute_cibersort`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_cibersort.md)
for the CIBERSORT R implementation
[`deconvolute_epic`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_epic.md)
for the EPIC algorithm
[`deconvolute_bayesdebulk`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_bayesdebulk.md)
for the BayesDeBulk algorithm

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
mixed_samples <- readRDS(data_file)

signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
signature_matrix <- readRDS(signature_file)

# Using EPIC algorithm
epic_results <- deconvolute(
  algorithm = "epic",
  data = mixed_samples,
  signature = signature_matrix,
  with_other_cells = TRUE
)

# Using BayesDeBulk algorithm
bayes_results <- deconvolute(
  algorithm = "bayesdebulk",
  data = mixed_samples,
  signature = signature_matrix,
  n_iter = 1000,
  burn_in = 100
)
} # }
```
