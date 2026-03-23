# Handle missing values

Imputes missing values in an expression matrix using various methods,
including simple replacement with the lowest non-zero value or more
sophisticated imputation techniques via the MsCoreUtils package.

## Usage

``` r
handle_missing_values(data, imputation_mode = "lowest_value", ...)
```

## Arguments

- data:

  A numeric matrix containing the data to be processed, with identifiers
  as row names and samples as columns.

- imputation_mode:

  A string specifying the imputation method to use:

  - `"lowest_value"` (default): Replaces NAs with the lowest non-zero
    value in the matrix

  - `"knn"`: k-nearest neighbor imputation

  - `"zero"`: Replace NAs with zeros

  - `"MLE"`: Maximum likelihood estimation

  - `"bpca"`: Bayesian principal component analysis

  - `"RF"`: Random Forest imputation

  - `"min"`: Replace NAs with minimum value in each column

  - `"MinDet"`: Deterministic minimum value imputation

  - `"MinProb"`: Probabilistic minimum value imputation

  - `"QRILC"`: Quantile regression imputation of left-censored data

  - `"mixed"`: Mixed imputation based on feature-wise missingness

  - `"nbavg"`: Impute with average of neighbors

- ...:

  Additional arguments passed to
  [`MsCoreUtils::impute_matrix`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html).
  See the documentation of that function for method-specific parameters.

## Value

A matrix with the same dimensions as the input, but with missing values
imputed according to the specified method.

## See also

[`impute_matrix`](https://rdrr.io/pkg/MsCoreUtils/man/imputation.html)
for detailed description of the imputation methods

## Examples

``` r
# Create example matrix with missing values
mat <- matrix(c(1.2, 3.4, NA, 5.6, NA, 7.8, 9.0, 2.1, 4.3), nrow = 3, ncol = 3)
rownames(mat) <- c("Protein1", "Protein2", "Protein3")
colnames(mat) <- c("Sample1", "Sample2", "Sample3")

# View original matrix
print(mat)
#>          Sample1 Sample2 Sample3
#> Protein1     1.2     5.6     9.0
#> Protein2     3.4      NA     2.1
#> Protein3      NA     7.8     4.3

# Impute missing values with lowest non-zero value (default)
result1 <- handle_missing_values(mat)
print(result1)
#>          Sample1 Sample2 Sample3
#> Protein1     1.2     5.6     9.0
#> Protein2     3.4     1.2     2.1
#> Protein3     1.2     7.8     4.3

# Impute missing values using k-nearest neighbors
if (FALSE) { # \dontrun{
# Requires the 'impute' package
result2 <- handle_missing_values(mat, imputation_mode = "knn")
print(result2)
} # }
```
