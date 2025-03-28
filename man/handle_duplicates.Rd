% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_duplicates.R
\name{handle_duplicates}
\alias{handle_duplicates}
\title{Handle duplicate identifiers in an expression matrix}
\usage{
handle_duplicates(data, duplicate_mode = "slice")
}
\arguments{
\item{data}{A numeric matrix containing expression data with identifiers as
row names and samples as columns.}

\item{duplicate_mode}{A string specifying the approach to handle duplicates:
\itemize{
\item \code{"slice"}: Keep only the row with the maximum median value for each identifier (default)
\item \code{"merge"}: Merge duplicate rows by taking the column-wise median of values
}}
}
\value{
A numeric matrix with unique identifiers as row names. The number of
rows will be equal to the number of unique identifiers in the input matrix.
}
\description{
Resolves duplicate row identifiers in an expression matrix using the
specified method.
}
\examples{
# Create example matrix with duplicate identifiers
mat <- matrix(1:12, nrow = 4, ncol = 3)
rownames(mat) <- c("ID1", "ID2", "ID1", "ID3")
colnames(mat) <- c("Sample1", "Sample2", "Sample3")

# View original matrix
print(mat)

# Handle duplicates by keeping rows with maximum median (default)
result1 <- handle_duplicates(mat)
print(result1)

# Handle duplicates by merging rows
result2 <- handle_duplicates(mat, duplicate_mode = "merge")
print(result2)

}
