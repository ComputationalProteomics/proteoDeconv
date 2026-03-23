# Create signature matrix using CIBERSORTx

Generates a signature matrix from reference profiles and cell type
classes using the CIBERSORTx Docker image. This signature matrix
contains cell type-specific marker proteins that will be used for
deconvolution.

## Usage

``` r
create_signature_matrix(
  refsample = NULL,
  phenoclasses = NULL,
  g_min = 200,
  g_max = 400,
  q_value = 0.01,
  replicates = 5,
  sampling = 0.5,
  fraction = 0.75,
  filter = FALSE,
  verbose = FALSE,
  QN = FALSE,
  single_cell = FALSE,
  use_sudo = FALSE,
  ...
)
```

## Arguments

- refsample:

  A numeric matrix containing reference profiles with genes as row names
  and samples as columns. This should be preprocessed data from pure
  cell populations.

- phenoclasses:

  A numeric matrix, data frame, or tibble containing cell type
  classification (0/1/2). If provided as a matrix or base data frame,
  the cell types are assumed to be in the row names. If provided as a
  tibble, the cell type identifiers should be included as a column. This
  is typically created using the create_phenoclasses() function.

- g_min:

  Minimum number of genes per cell type in the signature matrix. Default
  is 200.

- g_max:

  Maximum number of genes per cell type in the signature matrix. Default
  is 400.

- q_value:

  Q-value threshold for differential expression. Default is 0.01.

- replicates:

  Number of replicates to use for building the reference file (only
  relevant when single_cell=TRUE). Default is 5.

- sampling:

  Fraction of available gene expression profiles selected by random
  sampling (only relevant when single_cell=TRUE). Default is 0.5.

- fraction:

  Fraction of cells of the same identity showing evidence of expression
  (only relevant when single_cell=TRUE). Default is 0.75.

- filter:

  Logical indicating whether to remove non-hematopoietic genes. Default
  is FALSE.

- verbose:

  Logical indicating whether to print detailed output. Default is FALSE.

- QN:

  Logical indicating whether to run quantile normalization. Default is
  FALSE.

- single_cell:

  Logical indicating whether to create signature from scRNA-Seq data.
  Default is FALSE.

- use_sudo:

  Logical indicating whether to use sudo for Docker commands. Default is
  FALSE.

- ...:

  Additional arguments passed to the function.

## Value

A numeric matrix with genes as rows and cell types as columns,
representing the expression profile of each cell type. This signature
matrix can be used as input for deconvolution algorithms to estimate
cell type proportions in mixed samples.

## Details

This function uses the CIBERSORTx Docker image to construct a signature
matrix based on the input reference profiles and cell type
classifications. The CIBERSORTx token and email must be set as
environment variables either in the project directory's .Renviron file
or in the user's home directory .Renviron file.

## See also

[`create_phenoclasses`](https://computationalproteomics.github.io/proteoDeconv/reference/create_phenoclasses.md)
for creating the phenoclasses input and
[`deconvolute`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute.md)
for using the signature matrix in deconvolution.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load preprocessed pure samples data
pure_samples <- readRDS(system.file("extdata", "pure_samples_matrix.rds",
                                    package = "proteoDeconv"))

# Create phenoclasses for the samples
mapping_rules <- list(
  "CD8+ T cells" = "CD8",
  "Monocytes" = "Mono"
)
phenoclasses <- create_phenoclasses(
  data = pure_samples,
  mapping_rules = mapping_rules
)

# Create signature matrix
signature_matrix <- create_signature_matrix(
  refsample = pure_samples,
  phenoclasses = phenoclasses,
  g_min = 200,
  g_max = 400,
  q_value = 0.01
)
} # }
```
