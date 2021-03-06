% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complex_functions.R
\name{test_gene_list_in_exposures}
\alias{test_gene_list_in_exposures}
\title{Test if mutated PIDs are enriched in signatures}
\usage{
test_gene_list_in_exposures(
  in_gene_list,
  in_exposure_df,
  in_mut_table,
  in_gene.field = "GENE_short",
  in_p_cutoff = 0.05
)
}
\arguments{
\item{in_gene_list}{List with genes of interest}

\item{in_exposure_df}{Data frame with the signature exposures}

\item{in_mut_table}{Data frame or table of mutations (derived from
vcf-format)}

\item{in_gene.field}{Name of the column in which the gene names are to be
looked up}

\item{in_p_cutoff}{Significance threshold}
}
\value{
A list with entries \code{pvals}, \code{exposure_df},
  \code{number_of_mutated}, \itemize{ \item \code{pvals}: p-values of the
  t-tests performed on mutated vs. unmutated PIDs \item \code{exposure_df}:
  Transposed input exposures data frame with additional annotations for
  mutation status \item \code{number_of_mutated}: Number of PIDs carrying a
  mutation }
}
\description{
For all signatures found in a project, this function tests whether PIDs
having mutations in a specified list of genes of interest have significantly
higher exposures.
}
\examples{
 NULL

}
