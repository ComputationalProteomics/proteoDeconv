% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolute_epic.R
\name{deconvolute_epic}
\alias{deconvolute_epic}
\title{Deconvolute Bulk Data using EPIC}
\usage{
deconvolute_epic(
  preprocessed_data,
  signature_df,
  with_other_cells = TRUE,
  method_label
)
}
\arguments{
\item{preprocessed_data}{A matrix or data frame of the bulk proteome data.}

\item{signature_df}{A matrix or data frame of reference signature profiles.}

\item{with_other_cells}{Logical; if TRUE EPIC will include other cell types. Default is TRUE.}

\item{method_label}{A character string label used to tag the output with the deconvolution method.}
}
\value{
A tibble in long format with columns for \code{sample}, \code{cell_type}, \code{cell_count}, and \code{method} (with additional fields as provided in \code{...}).
}
\description{
Runs the EPIC algorithm on preprocessed proteomic data to estimate cell type fractions.
}
\details{
The function preprocesses both the bulk data and the signature matrix using \code{handle_scaling}
and \code{handle_input_data}. It then runs the EPIC deconvolution method from the EPIC package and reshapes
the resulting mRNA proportions into a tidy format.
}
\examples{
\dontrun{
result <- deconvolute_epic(bulk_data, signature_matrix, with_other_cells = TRUE, method_label = "EPIC")
}

}
