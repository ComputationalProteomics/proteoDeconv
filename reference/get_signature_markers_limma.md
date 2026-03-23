# Get signature markers using limma-based fold change method

Extracts unique signature markers from a signature matrix using a
fold-change approach. This method calculates fold changes between cell
types and selects genes that have significant differential expression.

## Usage

``` r
get_signature_markers_limma(
  Y,
  signature_mat,
  min_fold_change = 3,
  min_expression = 0,
  pseudo_count = 1e-10
)
```

## Arguments

- Y:

  A numeric matrix of bulk data with gene identifiers as row names.

- signature_mat:

  A numeric matrix of signature marker values with gene identifiers as
  row names and cell types as columns.

- min_fold_change:

  Minimum fold change required for a gene to be considered a marker
  (default: 3).

- min_expression:

  Minimum expression level required in the target cell type (default:
  0).

- pseudo_count:

  Small value added to expression values to avoid division by zero
  (default: 1e-10).

## Value

A character matrix of markers with columns for cell type, comparison
cell type, and gene name.
