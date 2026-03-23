# Get signature markers using simple threshold method

Extracts unique signature markers from a signature matrix using a simple
threshold-based approach. A gene is considered a marker if its
expression in one cell type exceeds its expression in other cell types
by a specified fold change and meets a minimum expression threshold.

## Usage

``` r
get_signature_markers_simple(
  Y,
  signature_mat,
  min_expression = 1000,
  min_fold_change = 3
)
```

## Arguments

- Y:

  A numeric matrix of bulk data with gene identifiers as row names.

- signature_mat:

  A numeric matrix of signature marker values with gene identifiers as
  row names and cell types as columns.

- min_expression:

  Minimum expression level required in the target cell type (default:
  1000).

- min_fold_change:

  Minimum fold change required compared to other cell types (default:
  3).

## Value

A character matrix of markers with columns for cell type, comparison
cell type, and gene name.
