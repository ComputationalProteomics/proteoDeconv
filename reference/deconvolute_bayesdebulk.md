# Deconvolute bulk proteome data using BayesDeBulk

Deconvolutes bulk proteome data using the [BayesDeBulk
algorithm](https://github.com/WangLab-MSSM/BayesDeBulk) to estimate cell
type proportions in mixed samples based on a signature matrix.

## Usage

``` r
deconvolute_bayesdebulk(
  data,
  signature,
  n_iter = 1000,
  burn_in = 100,
  marker_selection = "limma",
  ...
)
```

## Arguments

- data:

  A numeric matrix of bulk proteome data with gene identifiers as row
  names and samples as columns.

- signature:

  A numeric matrix containing signature marker values with gene
  identifiers as row names and cell types as columns.

- n_iter:

  Number of iterations for the MCMC sampling; default is 1000.

- burn_in:

  Number of burn-in iterations to discard; default is 100.

- marker_selection:

  The method to use for marker selection: "limma" (default), "simple",
  or a pre-computed marker matrix.

- ...:

  Additional arguments passed to marker selection functions.

## Value

A matrix containing cell type proportions with samples as rows and cell
types as columns.

## Details

This function calculates signature markers using the specified marker
selection method and runs the [BayesDeBulk deconvolution
algorithm](https://github.com/WangLab-MSSM/BayesDeBulk). The marker
selection process identifies genes that are uniquely expressed in
specific cell types, which are then used for the deconvolution.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data and signature matrix
data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
mixed_samples <- readRDS(data_file)

signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
signature_matrix <- readRDS(signature_file)

# Run deconvolution
result <- deconvolute_bayesdebulk(mixed_samples, signature_matrix)
} # }
```
