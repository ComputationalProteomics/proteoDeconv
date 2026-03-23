# Handle duplicate identifiers in an expression matrix

Resolves duplicate row identifiers in an expression matrix using the
specified method.

## Usage

``` r
handle_duplicates(data, duplicate_mode = "slice")
```

## Arguments

- data:

  A numeric matrix containing expression data with identifiers as row
  names and samples as columns.

- duplicate_mode:

  A string specifying the approach to handle duplicates:

  - `"slice"`: Keep only the row with the maximum median value for each
    identifier (default)

  - `"merge"`: Merge duplicate rows by taking the column-wise median of
    values

## Value

A numeric matrix with unique identifiers as row names. The number of
rows will be equal to the number of unique identifiers in the input
matrix.

## Examples

``` r
# Create example matrix with duplicate identifiers
mat <- matrix(1:12, nrow = 4, ncol = 3)
rownames(mat) <- c("ID1", "ID2", "ID1", "ID3")
colnames(mat) <- c("Sample1", "Sample2", "Sample3")

# View original matrix
print(mat)
#>     Sample1 Sample2 Sample3
#> ID1       1       5       9
#> ID2       2       6      10
#> ID1       3       7      11
#> ID3       4       8      12

# Handle duplicates by keeping rows with maximum median (default)
result1 <- handle_duplicates(mat)
print(result1)
#>     Sample1 Sample2 Sample3
#> ID1       3       7      11
#> ID2       2       6      10
#> ID3       4       8      12

# Handle duplicates by merging rows
result2 <- handle_duplicates(mat, duplicate_mode = "merge")
print(result2)
#>     Sample1 Sample2 Sample3
#> ID1       2       6      10
#> ID2       2       6      10
#> ID3       4       8      12
```
