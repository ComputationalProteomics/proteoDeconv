
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteoDeconv

<!-- badges: start -->

[![R-CMD-check](https://github.com/ComputationalProteomics/proteoDeconv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ComputationalProteomics/proteoDeconv/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

proteoDeconv is an R package that facilitates deconvolution of bulk
sample proteomic data to estimate proportions of cell types.

## Installation

You can install the development version of proteoDeconv from GitHub
with:

``` r
# install.packages("devtools")
devtools::install_github("ComputationalProteomics/proteoDeconv")
```

## Usage Example

The following example demonstrates a basic workflow using proteoDeconv.
The input file report.pg_matrix.tsv in this case is from DIA-NN, and the
LM7c.txt signature matrix was downloaded from the [Decomprolute
repository](https://github.com/PNNL-CompBio/decomprolute/tree/main/signature_matrices).

``` r
library(proteoDeconv)
pg <-
  read_tsv("report.pg_matrix.tsv") |>
  update_gene_symbols() |>
  handle_gene_groups() |>
  handle_missing_values(imputation_mode = "lowest_value") |>
  handle_duplicate_genes(duplicate_mode = "slice") |>
  handle_scaling(unlog = FALSE, tpm = TRUE)
  
signature <- read_tsv("LM7c.txt")
pg_deconvoluted <- deconvolute("cibersortx", pg, signature)
```

## Prerequisites

To use proteoDeconv with CIBERSORTx, you will need to register on the
[CIBERSORTx website](https://cibersortx.stanford.edu) and request a
token. This token is required for accessing the CIBERSORTx functionality
and is free for academic users. The token and email used can be stored
as environment variables in ~/.Renviron:

    CIBERSORTX_TOKEN = your_token_here
    CIBERSORTX_EMAIL = your_email_here

Additionally, Docker needs to be installed on the system in order to run
CIBERSORTx.

For running proteoDeconv with the regular CIBERSORT (not CIBERSORTx),
you will need to download the CIBERSORT source code from the [CIBERSORTx
website](https://cibersortx.stanford.edu).

Additionally, you can retrieve ready-made signature matrices for cell
type deconvolution from the [Decomprolute repository on
GitHub](https://github.com/PNNL-CompBio/decomprolute/tree/main/signature_matrices).
These matrices can be used as inputs in the deconvolution process with
proteoDeconv.
