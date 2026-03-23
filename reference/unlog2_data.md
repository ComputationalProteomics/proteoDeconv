# Convert data from log2 scale to linear scale

Transforms log2-transformed expression data back to linear scale by
calculating 2^x for each value in the input matrix.

## Usage

``` r
unlog2_data(data)
```

## Arguments

- data:

  A numeric matrix containing log2-transformed data with identifiers as
  row names and samples as columns.

## Value

A numeric matrix with the same dimensions as the input, with values
converted to linear scale (2^x).

## Details

This function can be used to convert log2-transformed proteomics or gene
expression data back to their original scale for downstream analyses
that require linear values. It preserves the row and column names of the
input matrix.

## Examples

``` r
# Create log2-transformed data matrix
log2_mat <- matrix(rnorm(12, mean = 10, sd = 2), nrow = 4, ncol = 3)
rownames(log2_mat) <- paste0("Protein", 1:4)
colnames(log2_mat) <- paste0("Sample", 1:3)

# View original log2 values
print(log2_mat)
#>            Sample1  Sample2  Sample3
#> Protein1  6.780849 12.25661 11.75693
#> Protein2  7.900580 10.87000 10.70160
#> Protein3 14.104068 11.09767 10.09976
#> Protein4 10.352873 11.29484 11.67150

# Convert to linear scale
linear_mat <- unlog2_data(log2_mat)
print(linear_mat)
#>             Sample1  Sample2  Sample3
#> Protein1   109.9611 4893.377 3460.889
#> Protein2   238.9525 1871.531 1665.334
#> Protein3 17609.5198 2191.445 1097.313
#> Protein4  1307.7516 2512.376 3261.905
```
