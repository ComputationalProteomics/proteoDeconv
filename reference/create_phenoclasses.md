# Create phenoclasses matrix for cell type deconvolution

Generates a phenotype classification matrix from sample data for use in
cell type deconvolution workflows. This matrix maps samples to their
respective cell types.

## Usage

``` r
create_phenoclasses(
  data,
  mapping_rules,
  verbose = FALSE,
  return_format = c("tibble", "matrix")
)
```

## Arguments

- data:

  A numeric matrix of pure cell type profiles with genes as row names
  and samples as columns. Column names should contain identifiers that
  can be mapped to cell types.

- mapping_rules:

  Either: 1. A named list where names are cell type labels and values
  are regular expression patterns used to match sample column names
  (e.g., list("CD8+ T cells" = "CD8", "Monocytes" = "Mono")), OR 2. A
  character vector with the same length and order as colnames(data),
  directly specifying the cell type for each sample.

- verbose:

  Logical. If TRUE, displays additional messages during processing.

- return_format:

  A string specifying the return format: "matrix" or "tibble" (default).

## Value

If return_format is "matrix", a numeric matrix with cell type groups as
rows and samples as columns. If return_format is "tibble", a tibble with
a "cell_type" column and columns for each sample.

## Details

The function is particularly useful for preparing reference data for
signature matrix generation using CIBERSORTx. Each cell in the output is
encoded as:

- 0 for 'Unknown' mappings

- 1 if the sample belongs to the cell type in that row

- 2 if the sample belongs to a different cell type

## See also

[`create_signature_matrix`](https://computationalproteomics.github.io/proteoDeconv/reference/create_signature_matrix.md)
which uses the phenoclasses matrix to generate a cell type signature
matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example using a named list of regex patterns
pure_samples <- readRDS(system.file("extdata", "pure_samples_matrix.rds",
                                     package = "proteoDeconv"))

mapping_rules <- list(
  "CD8+ T cells" = "CD8",
  "Monocytes" = "Mono"
)

phenoclasses1 <- create_phenoclasses(
  data = pure_samples,
  mapping_rules = mapping_rules,
  verbose = TRUE
)

# Example using a character vector of direct cell type assignments
cell_types <- c("CD8+ T cells", "CD8+ T cells", "Monocytes", "Monocytes", "Unknown")
phenoclasses2 <- create_phenoclasses(
  data = pure_samples,
  mapping_rules = cell_types,
  verbose = TRUE
)
} # }
```
