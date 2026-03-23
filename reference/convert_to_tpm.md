# Convert protein abundance data to TPM-like normalization

Normalizes protein abundance data using a TPM-like approach (Transcripts
Per Million), adapting this RNA-seq normalization method for use with
proteomics data.

## Usage

``` r
convert_to_tpm(data)
```

## Arguments

- data:

  A numeric matrix containing protein abundance data with identifiers as
  row names and samples as columns. Values should be in linear scale,
  not log-transformed.

## Value

A numeric matrix with the same dimensions as the input, with values
normalized to a TPM-like scale (sum of each column equals 1 million).

## Details

This function applies a TPM-like normalization to proteomics data, where
each protein abundance value is scaled by the total abundance in the
sample and multiplied by 1 million.

If your input data is log-transformed, use
[`unlog2_data`](https://computationalproteomics.github.io/proteoDeconv/reference/unlog2_data.md)
first to convert it to linear scale before applying this normalization.

## Examples

``` r
# Create example protein abundance data matrix
prot_mat <- matrix(abs(rnorm(12, mean = 500, sd = 200)), nrow = 4, ncol = 3)
rownames(prot_mat) <- paste0("Protein", 1:4)
colnames(prot_mat) <- paste0("Sample", 1:3)

# View original values and column sums
print(prot_mat)
#>            Sample1  Sample2  Sample3
#> Protein1 219.99130 624.3105 451.1601
#> Protein2 551.06341 729.6823 443.4589
#> Protein3  12.54728 135.6365 389.2601
#> Protein4 498.88574 450.5349 625.7964
print(colSums(prot_mat))
#>  Sample1  Sample2  Sample3 
#> 1282.488 1940.164 1909.676 

# Convert to TPM-like normalization
tpm_mat <- convert_to_tpm(prot_mat)

# Verify that column sums equal 1 million
print(tpm_mat)
#>             Sample1   Sample2  Sample3
#> Protein1 171534.816 321782.31 236249.6
#> Protein2 429683.184 376093.06 232216.9
#> Protein3   9783.546  69909.79 203835.7
#> Protein4 388998.453 232214.84 327697.8
print(colSums(tpm_mat))
#> Sample1 Sample2 Sample3 
#>   1e+06   1e+06   1e+06 
```
