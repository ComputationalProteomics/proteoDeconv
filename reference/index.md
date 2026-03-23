# Package index

## Data preprocessing

Functions for handling and processing protein data.

- [`handle_duplicates()`](https://computationalproteomics.github.io/proteoDeconv/reference/handle_duplicates.md)
  : Handle duplicate identifiers in an expression matrix
- [`extract_identifiers()`](https://computationalproteomics.github.io/proteoDeconv/reference/extract_identifiers.md)
  : Extract primary identifiers from e.g. protein groups
- [`update_gene_symbols()`](https://computationalproteomics.github.io/proteoDeconv/reference/update_gene_symbols.md)
  : Update gene symbols in a matrix to approved HGNC nomenclature
- [`unlog2_data()`](https://computationalproteomics.github.io/proteoDeconv/reference/unlog2_data.md)
  : Convert data from log2 scale to linear scale
- [`convert_to_tpm()`](https://computationalproteomics.github.io/proteoDeconv/reference/convert_to_tpm.md)
  : Convert protein abundance data to TPM-like normalization
- [`handle_scaling()`](https://computationalproteomics.github.io/proteoDeconv/reference/handle_scaling.md)
  : Handle scaling of protein abundance data
- [`handle_missing_values()`](https://computationalproteomics.github.io/proteoDeconv/reference/handle_missing_values.md)
  : Handle missing values

## Signature matrix creation

Functions for working with signature matrix creation.

- [`map_cell_groups()`](https://computationalproteomics.github.io/proteoDeconv/reference/map_cell_groups.md)
  : Map column names to cell type groups using patterns
- [`create_phenoclasses()`](https://computationalproteomics.github.io/proteoDeconv/reference/create_phenoclasses.md)
  : Create phenoclasses matrix for cell type deconvolution
- [`create_signature_matrix()`](https://computationalproteomics.github.io/proteoDeconv/reference/create_signature_matrix.md)
  : Create signature matrix using CIBERSORTx

## Deconvolution

Functions for deconvoluting bulk data into cell type proportions.

- [`deconvolute()`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute.md)
  : Deconvolute protein data using multiple algorithms
- [`deconvolute_bayesdebulk()`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_bayesdebulk.md)
  : Deconvolute bulk proteome data using BayesDeBulk
- [`deconvolute_cibersort()`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_cibersort.md)
  : Deconvolute mixture data using CIBERSORT
- [`deconvolute_cibersortx()`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_cibersortx.md)
  : Deconvolute using CIBERSORTx Docker image
- [`deconvolute_epic()`](https://computationalproteomics.github.io/proteoDeconv/reference/deconvolute_epic.md)
  : Deconvolute bulk data using EPIC

## Simulation

Function for simulating artifical mixtures.

- [`simulate_data()`](https://computationalproteomics.github.io/proteoDeconv/reference/simulate_data.md)
  : Simulate artifical mixtures from bulk proteome measurements
