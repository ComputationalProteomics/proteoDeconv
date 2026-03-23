# Update gene symbols in a matrix to approved HGNC nomenclature

This function checks and updates gene symbols in a matrix's row names to
ensure they conform to the current approved HGNC (HUGO Gene Nomenclature
Committee) standards. Non-standard or outdated gene symbols are replaced
with their current official symbols.

## Usage

``` r
update_gene_symbols(data, verbose = FALSE)
```

## Arguments

- data:

  A numeric matrix containing gene expression or protein abundance data
  with gene identifiers as row names.

- verbose:

  A logical value indicating whether to print detailed messages about
  the number of approved, non-approved, and unmappable gene symbols.
  Default is FALSE.

## Value

A matrix with updated gene symbols as row names. Rows with unmappable
symbols (where no suggested symbol could be found) are removed from the
output.

## Details

The function uses the `HGNChelper` package to check and standardize gene
symbols based on a predefined gene symbols map. This is important for
ensuring consistency in gene naming across datasets and avoiding issues
with outdated or non-standard gene symbols.
