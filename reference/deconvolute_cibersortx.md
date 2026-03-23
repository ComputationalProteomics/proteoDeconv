# Deconvolute using CIBERSORTx Docker image

Performs deconvolution of bulk proteome data into constituent cell types
using the CIBERSORTx Docker image. This function handles the interaction
with the Docker container and processes the results.

## Usage

``` r
deconvolute_cibersortx(
  data,
  signature,
  perm = 1,
  rmbatch_S_mode = FALSE,
  source_GEPs = NULL,
  rmbatch_B_mode = FALSE,
  QN = FALSE,
  absolute = FALSE,
  abs_method = "sig.score",
  use_sudo = FALSE
)
```

## Arguments

- data:

  A numeric matrix containing mixture data with genes as row names and
  samples as columns.

- signature:

  A numeric matrix containing the signature matrix with genes as row
  names and cell types as columns.

- perm:

  An integer specifying the number of permutations to be performed.
  Default is 1.

- rmbatch_S_mode:

  A logical value indicating whether to remove batch effects in source
  GEPs mode. Default is FALSE.

- source_GEPs:

  A matrix containing the source gene expression profiles. Required if
  `rmbatch_S_mode` is TRUE.

- rmbatch_B_mode:

  A logical value indicating whether to remove batch effects in bulk
  mode. Default is FALSE.

- QN:

  A logical value indicating whether to perform quantile normalization.
  Default is FALSE.

- absolute:

  A logical value indicating whether to use absolute mode. Default is
  FALSE.

- abs_method:

  A character string specifying the method to use for absolute mode.
  Default is "sig.score".

- use_sudo:

  A logical value indicating whether to use sudo for Docker commands.
  Default is FALSE.

## Value

A numeric matrix with samples as rows and cell types as columns,
representing the estimated proportion of each cell type in each sample.

## Details

This function requires the CIBERSORTx Docker image to be installed and
the `CIBERSORTX_EMAIL` and `CIBERSORTX_TOKEN` environment variables to
be set. You can get these credentials by registering at the CIBERSORTx
website (https://cibersortx.stanford.edu/).

The function creates temporary files for the mixture data and signature
matrix, runs the CIBERSORTx Docker container, and processes the results.
Note that absolute mode is not currently supported in the Docker
version.

## See also

[`deconvolute_cibersort`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_cibersort.md)
for using the R implementation of CIBERSORT,
[`deconvolute`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute.md)
for a unified interface to multiple deconvolution methods.

## Examples

``` r
if (FALSE) { # \dontrun{
# Set required environment variables (ideally in .Renviron)
Sys.setenv(CIBERSORTX_EMAIL = "your.email@example.com")
Sys.setenv(CIBERSORTX_TOKEN = "your-token-here")

# Load example data and signature matrix
data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
mixed_samples <- readRDS(data_file)

signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
signature_matrix <- readRDS(signature_file)

# Run deconvolution with CIBERSORTx Docker
result <- deconvolute_cibersortx(
  data = mixed_samples,
  signature = signature_matrix,
  perm = 100,
  QN = TRUE
)
} # }
```
