% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_functions.R
\name{cosineDist}
\alias{cosineDist}
\title{Compute the cosine distance of two vectors}
\usage{
cosineDist(a, b)
}
\arguments{
\item{a, b}{Numerical vectors of same length}
}
\value{
The scalar product of the two input vectors divided by the
       product of the norms of the two input vectors
}
\description{
Compute the cosine distance of two vectors
}
\examples{
## 1. Orthogonal vectors:
cosineDist(c(1,0),c(0,1))
## 2. Non-orthogonal vectors:
cosineDist(c(1,0),c(1,1))
## Compare trigonometry:
1-cos(pi/4)

}
