% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_missing_values.R
\name{handle_missing_values}
\alias{handle_missing_values}
\title{Handle Missing Values in a Data Frame}
\usage{
handle_missing_values(
  data,
  gene_column = "Genes",
  imputation_mode = "lowest_value",
  ...
)
}
\arguments{
\item{data}{A data frame containing the data to be processed.}

\item{gene_column}{A string specifying the name of the column containing gene identifiers. Default is "Genes".}

\item{imputation_mode}{A string specifying the imputation method to use. Options are "lowest_value", "knn", "zero", "MLE", "bpca", "RF", "min", "MinDet", "MinProb", "QRILC", "mixed", "nbavg", or "with". Default is "lowest_value".}

\item{...}{Additional arguments passed to MsCoreUtils::impute_matrix.}
}
\value{
A data frame with missing values imputed.
}
\description{
This function handles missing values in a data frame by either replacing them with the lowest non-zero value
or using k-nearest neighbors (KNN) imputation.
}
