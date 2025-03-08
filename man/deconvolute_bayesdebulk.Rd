% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolute_bayesdebulk.R
\name{deconvolute_bayesdebulk}
\alias{deconvolute_bayesdebulk}
\title{Deconvolute bulk proteome data using BayesDeBulk}
\usage{
deconvolute_bayesdebulk(
  data,
  signature,
  n_iter = 1000,
  burn_in = 100,
  marker_selection = "limma",
  ...
)
}
\arguments{
\item{data}{A numeric matrix of bulk proteome data with gene identifiers as
row names and samples as columns.}

\item{signature}{A numeric matrix containing signature marker values with
gene identifiers as row names and cell types as columns.}

\item{n_iter}{Number of iterations for the MCMC sampling; default is 1000.}

\item{burn_in}{Number of burn-in iterations to discard; default is 100.}

\item{marker_selection}{The method to use for marker selection: "limma"
(default), "simple", or a pre-computed marker matrix.}

\item{...}{Additional arguments passed to marker selection functions.}
}
\value{
A matrix containing cell type proportions with samples as rows and
cell types as columns.
}
\description{
Deconvolutes bulk proteome data using the \href{https://github.com/WangLab-MSSM/BayesDeBulk}{BayesDeBulk algorithm} to estimate cell type
proportions in mixed samples based on a signature matrix.
}
\details{
This function calculates signature markers using the specified
marker selection method and runs the \href{https://github.com/WangLab-MSSM/BayesDeBulk}{BayesDeBulk deconvolution algorithm}. The marker
selection process identifies genes that are uniquely expressed in specific
cell types, which are then used for the deconvolution.
}
\examples{
\dontrun{
# Load example data and signature matrix
data_file <- system.file("extdata", "mixed_samples_matrix.rds", package = "proteoDeconv")
mixed_samples <- readRDS(data_file)

signature_file <- system.file("extdata", "cd8t_mono_signature_matrix.rds", package = "proteoDeconv")
signature_matrix <- readRDS(signature_file)

# Run deconvolution
result <- deconvolute_bayesdebulk(mixed_samples, signature_matrix)
}

}
