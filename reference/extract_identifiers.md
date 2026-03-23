# Extract primary identifiers from e.g. protein groups

Extracts the primary identifier from compound row identifiers in a
matrix by splitting on a separator and keeping only the first entry.

## Usage

``` r
extract_identifiers(data, separator = ";")
```

## Arguments

- data:

  A numeric matrix with compound identifiers (e.g. protein/gene groups)
  as row names.

- separator:

  A character string used to separate multiple identifiers. Default is
  ";".

## Value

A matrix with simplified identifiers as row names, where each identifier
is the first element from the original compound identifier.

## Examples

``` r
# Create matrix with compound identifiers (like protein groups)
mat <- matrix(1:9, nrow = 3, ncol = 3)
rownames(mat) <- c("P04637;P02340", "Q15796;O35182", "P01308;P01315;P01317")
colnames(mat) <- c("Sample1", "Sample2", "Sample3")

# View original matrix
print(mat)
#>                      Sample1 Sample2 Sample3
#> P04637;P02340              1       4       7
#> Q15796;O35182              2       5       8
#> P01308;P01315;P01317       3       6       9

# Extract primary identifiers
result <- extract_identifiers(mat)
print(result)
#>        Sample1 Sample2 Sample3
#> P04637       1       4       7
#> Q15796       2       5       8
#> P01308       3       6       9
```
