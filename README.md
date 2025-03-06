
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
coverage](https://codecov.io/gh/manszamore/proteoDeconv/graph/badge.svg)](https://app.codecov.io/gh/manszamore/proteoDeconv)
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

mix_data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
mixed_samples <- readRDS(mix_data_file)

# Preprocess the mix samples data using the package's pipeline
mixed_samples_preprocessed <- mixed_samples |>
  extract_identifiers() |>
  update_gene_symbols() |>
  handle_missing_values(imputation_mode = "lowest_value") |>
  handle_duplicates(duplicate_mode = "slice") |>
  convert_to_tpm()

# Load an example signature matrix
signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
signature_matrix <- readRDS(signature_file)

# Perform deconvolution
results <- deconvolute(
  algorithm = "epic",
  data = mixed_samples_preprocessed,
  signature = signature_matrix,
  with_other_cells = FALSE
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
