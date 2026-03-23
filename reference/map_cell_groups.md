# Map column names to cell type groups using patterns

Maps a character vector of column names to predefined cell type
categories by matching against patterns, with optional regex support.

## Usage

``` r
map_cell_groups(
  column_names,
  mapping_rules,
  default_group = "Unknown",
  verbose = FALSE,
  use_regex = TRUE
)
```

## Arguments

- column_names:

  A character vector of column names to be mapped to cell type groups.

- mapping_rules:

  A named list where names are the target cell type groups and values
  are character vectors of patterns that identify columns belonging to
  each group.

- default_group:

  A string that represents the default group name for columns not
  matching any pattern. Default is "Unknown".

- verbose:

  Logical. If TRUE, prints mapping process messages showing which
  columns were mapped to which groups and which columns remained
  unmapped. Default is FALSE.

- use_regex:

  Logical. If TRUE (default), treats patterns in mapping_rules as
  regular expressions. If FALSE, performs exact string matching.

## Value

A character vector with the same length as `column_names`, containing
the mapped cell type group for each column name. Columns that don't
match any pattern will be assigned the `default_group`.

## Details

This function iterates through the `mapping_rules` list and attempts to
match each column name against the patterns for each cell type group.
The first matching group in the order of the list will be assigned.

When `use_regex = TRUE`, patterns are treated as regular expressions and
matching is case-insensitive. When `use_regex = FALSE`, exact string
matching is performed, which is case-sensitive.

## Examples

``` r
# Example column names from a dataset
cols <- c("CD8_T_cell", "T_cell", "B_cell",
          "NK_cell", "Monocyte", "Unknown_cell")

# Define simple mapping rules for cell types
mapping <- list(
  "T_cell" = c("CD8_T_cell", "T_cell"),
  "B_cell" = c("B_cell"),
  "NK_cell" = c("NK_cell"),
  "Monocyte" = c("Monocyte")
)

# Map column names to cell types using exact matching
cell_types <- map_cell_groups(cols, mapping, use_regex = FALSE)
print(data.frame(column = cols, cell_type = cell_types))
#>         column cell_type
#> 1   CD8_T_cell    T_cell
#> 2       T_cell    T_cell
#> 3       B_cell    B_cell
#> 4      NK_cell   NK_cell
#> 5     Monocyte  Monocyte
#> 6 Unknown_cell   Unknown

# Define mapping rules using regex patterns
mapping_regex <- list(
  "T_cell" = c("CD8.*", "T_cell"),
  "B_cell" = c("B[-_]?cell"),
  "NK_cell" = c("NK[-_]?cell"),
  "Monocyte" = c("Mono.*")
)

# Map using regex patterns (default)
cell_types_regex <- map_cell_groups(cols, mapping_regex)
print(data.frame(column = cols, cell_type = cell_types_regex))
#>         column cell_type
#> 1   CD8_T_cell    T_cell
#> 2       T_cell    T_cell
#> 3       B_cell    B_cell
#> 4      NK_cell   NK_cell
#> 5     Monocyte  Monocyte
#> 6 Unknown_cell   Unknown
```
