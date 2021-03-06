% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMC.R
\name{SMC_perPID}
\alias{SMC_perPID}
\title{Run SMC at a per sample level}
\usage{
SMC_perPID(
  in_dfList,
  in_LCDlist,
  in_subgroups_df,
  in_save_plot = TRUE,
  in_save_dir = NULL,
  in_save_name = "KataegisSMCs.pdf",
  in_verbose_flag = 0,
  ...
)
}
\arguments{
\item{in_dfList}{Named list of vcf-like data frames, one entry per sample/PID
of a cohort.}

\item{in_LCDlist}{Output of an LCD list perfomed on the above cohort,
carrying notably information on the exposures
(\code{in_LCDlist$exposures}), the present signatures
(\code{in_LCDlist$signatures}) and meta information about the
signatures(\code{in_LCDlist$out_sig_ind_df}).}

\item{in_subgroups_df}{Data frame with subgroup information about the PIDs in
the above mentioned cohort.}

\item{in_save_plot}{Boolean flag to indicate whether per-PID plots should be
saved.}

\item{in_save_dir}{If per-PID plots are to be saved, this is the path where
to save them.}

\item{in_save_name}{Suffix to be appended to the sample name to generate the
name of the saved per-PID plots.}

\item{in_verbose_flag}{Whether to run verbose (1) or not (0).}

\item{...}{Data passed on to \code{\link{run_SMC}}.}
}
\value{
A list of lists. The top level is a named per-PID list, each entry is
  of type SMClist (cf. \code{\link{run_SMC}}).
}
\description{
Run an SMC analysis (stratification of the mutational catalogue) at per
sample / per-PID level, corresponding to a divide and conquer strategy. For
every single PID, only those signatures actually present in this PID will be
provided for the SMC analysis.
}
\examples{
 NULL
 
}
