% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_functions.R
\name{normalize_df_per_dim}
\alias{normalize_df_per_dim}
\alias{average_over_present}
\alias{sd_over_present}
\alias{stderrmean_over_present}
\title{Useful functions on data frames}
\usage{
normalize_df_per_dim(in_df, in_dimension)

average_over_present(in_df, in_dimension)

sd_over_present(in_df, in_dimension)

stderrmean_over_present(in_df, in_dimension)
}
\arguments{
\item{in_df}{Data frame to be normalized}

\item{in_dimension}{Dimension along which the operation will be carried out}
}
\value{
The normalized numerical data frame (\code{normalize_df_per_dim})

A vector of the means (\code{average_over_present})

A vector of the standard deviations (\code{sd_over_present})

A vector of the standard errors of the mean 
 (\code{stderrmean_over_present})
}
\description{
\code{normalize_df_per_dim}: Normalization is carried out by dividing by
\code{rowSums} or \code{colSums}; for rows with \code{rowSums=0} or columns
with \code{colSums=0}, the normalization is left out.

\code{average_over_present}: 
If averaging over columns, zero rows (i.e. those with \code{rowSums=0})
are left out, if averaging over rows, zero columns (i.e. those with
\code{colSums=0}) are left out.

\code{sd_over_present}: 
If computing the standard deviation over columns, zero rows
(i.e. those with \code{rowSums=0}) are left out, if computing
the standard deviation over rows, zero columns (i.e. those with
\code{colSums=0}) are left out.

\code{stderrmean_over_present}: 
If computing the standard error of the mean over columns, zero rows
(i.e. those with \code{rowSums=0}) are left out, if computing the
standard error of the mean over rows, zero columns (i.e. those with
\code{colSums=0}) are left out. Uses the function \code{\link{stderrmean}}
}
\examples{
test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),
                             ncol=4))
## 1. Normalize over rows:
normalize_df_per_dim(test_df,1)
## 2. Normalize over columns:
normalize_df_per_dim(test_df,2)

test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),
                             ncol=4))
## 1. Average over non-zero rows:
average_over_present(test_df,1)
## 2. Average over non-zero columns:
average_over_present(test_df,2)

test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),
                             ncol=4))
## 1. Compute standard deviation over non-zero rows:
sd_over_present(test_df,1)
## 2. Compute standard deviation over non-zero columns:
sd_over_present(test_df,2)

test_df <- data.frame(matrix(c(1,2,3,0,5,2,3,4,0,6,0,0,0,0,0,4,5,6,0,7),
                             ncol=4))
## 1. Compute standard deviation over non-zero rows:
stderrmean_over_present(test_df,1)
## 2. Compute standard deviation over non-zero columns:
stderrmean_over_present(test_df,2)

}
\seealso{
\code{\link{stderrmean}}
}
