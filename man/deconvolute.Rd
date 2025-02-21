% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolute.R
\name{deconvolute}
\alias{deconvolute}
\title{Deconvolute Protein Data}
\usage{
deconvolute(algorithm, data, signature_matrix = NULL, ...)
}
\arguments{
\item{algorithm}{A string specifying the deconvolution algorithm to use. Options are "cibersortx", "cibersort", "epic", and other algorithms supported by the \code{immunedeconv} package.}

\item{data}{A data frame containing protein or gene data. The data frame should have proteins/genes as rows and samples as columns.}

\item{signature_matrix}{An optional signature matrix required by some algorithms.}

\item{...}{Additional arguments passed to the specific deconvolution functions.}
}
\value{
A data frame with the deconvolution results.
}
\description{
This function deconvolutes protein or gene expression data using various algorithms.
}
\details{
The function supports multiple deconvolution algorithms:
\itemize{
\item \strong{"cibersortx"}: Uses the \code{deconvolute_cibersortx} function.
\item \strong{"cibersort"}: Uses the \code{deconvolute_cibersort} function.
\item \strong{"epic"}: Uses the \code{deconvolute_epic} function.
}
}
