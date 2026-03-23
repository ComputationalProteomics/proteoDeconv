# Get signature markers for deconvolution

Extracts unique signature markers from a signature matrix for use in
cell type deconvolution. This function selects a method to identify cell
type-specific markers.

## Usage

``` r
get_signature_markers(Y, signature_mat, method = "limma", ...)
```

## Arguments

- Y:

  A numeric matrix of bulk data with gene identifiers as row names.

- signature_mat:

  A numeric matrix of signature marker values with gene identifiers as
  row names and cell types as columns.

- method:

  The method to use for marker selection: "limma" (default), "simple",
  or a pre-computed marker matrix.

- ...:

  Additional arguments passed to the selected marker selection function.

## Value

A character matrix of markers with columns for cell type, comparison
cell type, and gene name.
