% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolute_bayesdebulk.R
\name{get_signature_markers_limma}
\alias{get_signature_markers_limma}
\title{Get signature markers using limma-based fold change method}
\usage{
get_signature_markers_limma(
  Y,
  signature_mat,
  min_fold_change = 3,
  min_expression = 0,
  pseudo_count = 1e-10
)
}
\arguments{
\item{Y}{A numeric matrix of bulk data with gene identifiers as row names.}

\item{signature_mat}{A numeric matrix of signature marker values with gene identifiers as row names
and cell types as columns.}

\item{min_fold_change}{Minimum fold change required for a gene to be considered a marker (default: 3).}

\item{min_expression}{Minimum expression level required in the target cell type (default: 0).}

\item{pseudo_count}{Small value added to expression values to avoid division by zero (default: 1e-10).}
}
\value{
A character matrix of markers with columns for cell type, comparison cell type, and gene name.
}
\description{
Extracts unique signature markers from a signature matrix using a fold-change approach.
This method calculates fold changes between cell types and selects genes that have
significant differential expression.
}
\keyword{internal}
