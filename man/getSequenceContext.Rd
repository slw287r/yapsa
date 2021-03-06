% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indel_functions.R
\name{getSequenceContext}
\alias{getSequenceContext}
\title{Extracts the sequence context up and downstream of a nucleotide position}
\usage{
getSequenceContext(position, chr, offsetL = 10, offsetR = 50)
}
\arguments{
\item{position}{Start position of the considered INDEL}

\item{chr}{Chromosome of the considered INDEL}

\item{offsetL}{Number of nucleotides downstream of \code{position}}

\item{offsetR}{Number of nucleotides upstream of \code{position}}
}
\value{
Returns a character string containing the defined seqeunce context
}
\description{
Extracts the sequence context up and downstream of a nucleotide position
}
\examples{
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

sequence_context <- getSequenceContext(position = 123456789, chr = "chr12",
                                       offsetL= 10, offsetR=50)
sequence_context

}
