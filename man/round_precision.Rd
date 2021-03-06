% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_functions.R
\name{round_precision}
\alias{round_precision}
\title{Round to a defined precision}
\usage{
round_precision(x, breaks = NULL, in_precision = 0.05, in_interval = c(0, 1))
}
\arguments{
\item{x}{Vector to be rounded}

\item{breaks}{The breaks used for rounding. Default NULL}

\item{in_precision}{Precition default 0.05}

\item{in_interval}{Interval needs to be larger than the precision value}
}
\value{
A list with two entries: \itemize{ \item \code{values}: the rounded
  vector \item \code{breaks}: the breaks used for rounding }
}
\description{
This function is an extension with regard to the function \code{\link{round}}
from base R as it allows not only digits as precision, but can also round to
a user-specified precision. The interval in which the rounding operation is
to be carried out also can be specified by the user (default is the unit
interval). Alternatively, breaks can be provided.
}
\examples{
 NULL

}
