% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMC.R
\name{run_SMC}
\alias{run_SMC}
\title{Wrapper function for the Stratification of a Mutational Catalogue}
\usage{
run_SMC(
  my_table,
  this_signatures_df,
  this_signatures_ind_df,
  this_subgroups_df,
  column_name,
  refGenome,
  cohort_method_flag = "all_PIDs",
  in_strata_order_ind = seq_len(length(unique(my_table[, column_name]))),
  wordLength = 3,
  verbose_flag = 1,
  target_dir = NULL,
  strata_dir = NULL,
  output_path = NULL,
  in_all_exposures_df = NULL,
  in_rownames = c(),
  in_norms = NULL,
  in_label_orientation = "turn",
  this_sum_ind = NULL
)
}
\arguments{
\item{my_table}{A big dataframe constructed from a vcf-like file of a whole
cohort. The first columns are those of a standard vcf file, followed by an
arbitrary number of custom or user defined columns. One of these must carry
a PID (patient or sample identifyier) and one must be the category used for
stratification.}

\item{this_signatures_df}{A numeric data frame \code{W} in with \code{n} rows
and \code{l} columns, \code{n} being the number of features and \code{l}
being the number of signatures}

\item{this_signatures_ind_df}{A data frame containing meta information about
the signatures}

\item{this_subgroups_df}{A data frame indicating which PID (patient or sample
identifyier) belongs to which subgroup}

\item{column_name}{Name of the column in \code{my_table} which is going to be
used for stratification}

\item{refGenome}{FaFile of the reference genome to extract the motif context
of the variants in \code{my_table}}

\item{cohort_method_flag}{Either or several of
\code{c("all_PIDs","cohort","norm_PIDs")}, representing alternative ways to
average over the cohort.}

\item{in_strata_order_ind}{Index vector defining reordering of the strata}

\item{wordLength}{Integer number defining the length of the features or
motifs, e.g. 3 for tripletts or 5 for pentamers}

\item{verbose_flag}{Verbose if \code{verbose_flag=1}}

\item{target_dir}{Path to directory where the results of the stratification
procedure are going to be stored if non-NULL.}

\item{strata_dir}{Path to directory where the mutational catalogues of the
different strata are going to be stored if non-NULL}

\item{output_path}{Path to directory where the results, especially the
figures produced by \code{\link{plot_SMC}} are going to be stored.}

\item{in_all_exposures_df}{Optional argument, if specified, \code{H}, i.e.
the overall exposures without stratification, is set to equal
\code{in_all_exposures_df}. This is equivalent to forcing the
\code{\link{LCD_SMC}} procedure to use e.g. the exposures of a previously
performed NMF decomposition.}

\item{in_rownames}{Optional parameter to specify rownames of the mutational
catalogue \code{V} i.e. the names of the features.}

\item{in_norms}{If specified, vector of the correction factors for every
motif due to differing trinucleotide content. If null, no correction is
applied.}

\item{in_label_orientation}{Whether or not to turn the labels on the x-axis.}

\item{this_sum_ind}{Optional set of indices for reordering the PIDs}
}
\value{
A list with entries \code{exposures_list}, \code{catalogues_list},
  \code{cohort} and \code{name_list}. \itemize{ \item \code{exposures_list}:
  The list of \code{s} strata specific exposures Hi, all are numerical data
  frames with \code{l} rows and \code{m} columns, \code{l} being the number
  of signatures and \code{m} being the number of samples \item
  \code{catalogues_list}: A list of \code{s} strata specific cohortwide (i.e.
  averaged over cohort) normalized exposures \item \code{cohort}:
  \code{subgroups_df} adjusted for plotting \item \code{name_list}: Names of
  the contructed strata. }
}
\description{
\code{\link{run_SMC}} takes as input a big dataframe constructed from a
vcf-like file of a whole cohort. This wrapper function calls custom functions
to construct a mutational catalogue and stratify it according to categories
indicated by a special column in the input dataframe: \itemize{ \item
\code{\link{create_mutation_catalogue_from_df}} \item
\code{adjust_number_of_columns_in_list_of_catalogues} } This stratification
yields a collection of stratified mutational catalogues, these are
reformatted and sent to the custom function \code{\link{SMC}} and thus
indirectly to \code{\link{LCD_SMC}} to perform a signature analysis of the
stratified mutational catalogues. The result is then handed over to
\code{\link{plot_SMC}} for visualization.
}
\examples{
 library(BSgenome.Hsapiens.UCSC.hg19)
 data(sigs)
 data(lymphoma_test)
 data(lymphoma_cohort_LCD_results)
 strata_list <-
   cut_breaks_as_intervals(lymphoma_test_df$random_norm,
                           in_outlier_cutoffs=c(-4,4),
                           in_cutoff_ranges_list=list(c(-2.5,-1.5),
                                                      c(0.5,1.5)),
                           in_labels=c("small","intermediate","big"))
 lymphoma_test_df$random_cat <- strata_list$category_vector
 choice_ind <- (names(lymphoma_Nature2013_COSMIC_cutoff_exposures_df)
                \%in\% unique(lymphoma_test_df$PID))
 lymphoma_test_exposures_df <-
   lymphoma_Nature2013_COSMIC_cutoff_exposures_df[,choice_ind]
 temp_subgroups_df <- make_subgroups_df(lymphoma_test_df,
                                        lymphoma_test_exposures_df)
 mut_density_list <- run_SMC(lymphoma_test_df,
                             AlexCosmicValid_sig_df,
                             AlexCosmicValid_sigInd_df,
                             temp_subgroups_df,
                             column_name="random_cat",
                             refGenome=BSgenome.Hsapiens.UCSC.hg19,
                             cohort_method_flag="norm_PIDs",
                             in_rownames = rownames(AlexCosmicValid_sig_df))

}
\seealso{
\code{\link{create_mutation_catalogue_from_df}}

\code{\link{normalizeMotifs_otherRownames}}

\code{\link{plot_SMC}}
}
