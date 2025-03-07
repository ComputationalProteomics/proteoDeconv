
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteoDeconv

proteoDeconv is an R package that facilitates the deconvolution of bulk
proteomic data to estimate cell type proportions. With proteoDeconv, you
can preprocess your data to prepare it for deconvolution, create cell
type signature matrices, and perform deconvolution using multiple
algorithms.

<!-- badges: start -->

[![R-CMD-check](https://github.com/ComputationalProteomics/proteoDeconv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ComputationalProteomics/proteoDeconv/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ComputationalProteomics/proteoDeconv/graph/badge.svg)](https://app.codecov.io/gh/ComputationalProteomics/proteoDeconv)
<!-- badges: end -->

## Installation

You can install the development version of proteoDeconv from GitHub:

``` r
# install.packages("pak")
pak::pak("ComputationalProteomics/proteoDeconv")
```

## Usage example

Below is a brief example that demonstrates loading the provided example
data and signature matrix from the package and performing deconvolution.
You can learn more in `vignette("proteoDeconv")`.

``` r
library(proteoDeconv)

# Load example data
mix_data <- readRDS(system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv"))
signature <- readRDS(system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv"))

# Preprocess data
mix_processed <- mix_data |>
  extract_identifiers() |>  # Extract IDs
  update_gene_symbols() |>  # Update to current gene symbols
  handle_missing_values() |> # Handle NAs
  handle_duplicates() |> # Handle duplicates
  convert_to_tpm()          # Scale to TPM-like scale

# Run deconvolution
results <- deconvolute(
  algorithm = "epic",
  data = mix_processed,
  signature = signature
)

# View the deconvolution results
print(results)
```

## Prerequisites

### CIBERSORTx

To use proteoDeconv with CIBERSORTx:

- Ensure Docker is installed on your system.

- Register and obtain a token from the [CIBERSORTx
  website](https://cibersortx.stanford.edu).

- Set the token and email as environment variables in your .Renviron
  file (this file can be in your home directory or in your project
  folder):

  ``` r
  CIBERSORTX_TOKEN=your_token_here
  CIBERSORTX_EMAIL=your_email_here
  ```

### CIBERSORT

For running proteoDeconv with the original CIBERSORT method, download
the CIBERSORT source code from the [CIBERSORT
website](https://cibersortx.stanford.edu) and source it before executing
deconvolutions:

``` r
source("/path/to/CIBERSORT.R")
```

## Further reference

For a detailed walkthrough of a complete workflow, see
`vignette("proteoDeconv")` or [the package
website](https://computationalproteomics.github.io/proteoDeconv/).

## Support

For bug reports or any other inquiries, please [open an
issue](https://github.com/ComputationalProteomics/proteoDeconv/issues)
on our GitHub repository.
