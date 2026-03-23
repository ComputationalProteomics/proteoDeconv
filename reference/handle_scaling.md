# Handle scaling of protein abundance data

Scales protein abundance data by optionally converting from log2 scale
to linear scale and/or applying TPM-like normalization for proteomics
data comparison.

## Usage

``` r
handle_scaling(data, unlog = FALSE, tpm = FALSE)
```

## Arguments

- data:

  A numeric matrix containing protein abundance data with identifiers as
  row names and samples as columns.

- unlog:

  A logical value indicating whether to unlog the data (convert from
  log2 scale). Default is FALSE.

- tpm:

  A logical value indicating whether to apply TPM-like normalization
  (adapting Transcripts Per Million for proteomics). Default is FALSE.

## Value

A numeric matrix with the scaled protein abundance data according to the
specified options.

## Details

This function combines the functionality of
[`unlog2_data`](https://computationalproteomics.github.io/proteoDeconv/reference/unlog2_data.md)
and
[`convert_to_tpm`](https://computationalproteomics.github.io/proteoDeconv/reference/convert_to_tpm.md)
to provide a flexible way to handle common scaling operations for
proteomics data.

## See also

[`unlog2_data`](https://computationalproteomics.github.io/proteoDeconv/reference/unlog2_data.md)
for just converting from log2 scale,
[`convert_to_tpm`](https://computationalproteomics.github.io/proteoDeconv/reference/convert_to_tpm.md)
for just applying TPM-like normalization

## Examples

``` r
# Create log2-transformed protein abundance matrix
log2_mat <- matrix(rnorm(12, mean = 10, sd = 2), nrow = 4, ncol = 3)
rownames(log2_mat) <- paste0("Protein", 1:4)
colnames(log2_mat) <- paste0("Sample", 1:3)

# Convert from log2 to linear scale only
linear_mat <- handle_scaling(log2_mat, unlog = TRUE, tpm = FALSE)

# Convert from log2 to linear scale and then apply TPM-like normalization
tpm_mat <- handle_scaling(log2_mat, unlog = TRUE, tpm = TRUE)

# Apply TPM-like normalization to already linear protein abundance data
linear_data <- matrix(abs(rnorm(12, mean = 500, sd = 200)), nrow = 4, ncol = 3)
tpm_only <- handle_scaling(linear_data, unlog = FALSE, tpm = TRUE)
```
