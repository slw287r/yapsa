# Copyright © 2014-2019  The YAPSA package contributors
# This file is part of the YAPSA package. The YAPSA package is licenced under
# GPL-3

#'Linear Combination Decomposition
#'
#'\code{LCD} performs a mutational signatures decomposition of a given
#'mutational catalogue \code{V} with known signatures \code{W} by solving the
#'minimization problem \eqn{min(||W*H - V||)} with additional constraints of
#'non-negativity on H where W and V are known
#'
#'@param in_mutation_catalogue_df A numeric data frame \code{V} with \code{n}
#'  rows and \code{m} columns, \code{n} being the number of features and
#'  \code{m} being the number of samples
#'@param in_signatures_df A numeric data frame \code{W} with \code{n} rows and
#'  \code{l} columns, \code{n} being the number of features and \code{l} being
#'  the number of signatures
#'@param in_per_sample_cutoff A numeric value less than 1. Signatures from
#'  within \code{W} with an exposure per sample less than \code{in_cutoff} will
#'  be discarded.
#'
#'@return The exposures \code{H}, a numeric data frame with \code{l} rows and
#'  \code{m} columns, \code{l} being the number of signatures and \code{m} being
#'  the number of samples
#'
#'@seealso \code{\link[lsei]{lsei}}
#'
#' @examples
#'
#' ## define raw data
#' W_prim <- matrix(c(1,2,3,4,5,6),ncol=2)
#' W_prim_df <- as.data.frame(W_prim)
#' W_df <- YAPSA:::normalize_df_per_dim(W_prim_df,2) # corresponds to the sigs
#' W <- as.matrix(W_df)
#' ## 1. Simple case: non-negativity already in raw data
#' H <- matrix(c(2,5,3,6,1,9,1,2),ncol=4)
#' H_df <- as.data.frame(H) # corresponds to the exposures
#' V <- W %*% H # matrix multiplication
#' V_df <- as.data.frame(V) # corresponds to the mutational catalogue
#' exposures_df <- YAPSA:::LCD(V_df,W_df)
#' ## 2. more complicated: raw data already contains negative elements
#' ## define indices where sign is going to be swapped
#' sign_ind <- c(5,7)
#' ## now compute the indices of the other fields in the columns affected
#' ## by the sign change
#' row_ind <- sign_ind %% dim(H)[1]
#' temp_ind <- 2*row_ind -1
#' other_ind <- sign_ind + temp_ind
#' ## alter the matrix H to yield a new mutational catalogue
#' H_compl <- H
#' H_compl[sign_ind] <- (-1)*H[sign_ind]
#' H_compl_df <- as.data.frame(H_compl) # corresponds to the exposures
#' V_compl <- W %*% H_compl # matrix multiplication
#' V_compl_df <- as.data.frame(V_compl) # corresponds to the mutational catalog
#' exposures_df <- YAPSA:::LCD(V_compl_df,W_df)
#' exposures <- as.matrix(exposures_df)
#'
#'@importFrom lsei lsei
#'
#'@export
#'
LCD <- function(in_mutation_catalogue_df,
                in_signatures_df,
                in_per_sample_cutoff=0){
  signatures_matrix <- as.matrix(in_signatures_df)
  out_exposures_df <- data.frame()
  G <- diag(dim(signatures_matrix)[2])
  H <- rep(0,dim(signatures_matrix)[2])
  for (i in seq_len(ncol(in_mutation_catalogue_df))) {
    # temp_fractions <- limSolve::lsei(A = signatures_matrix,
    #                                  B = in_mutation_catalogue_df[,i],
    #                                  G=G, H=H, verbose=FALSE)
    # temp_exposures_vector <- as.vector(temp_fractions$X)
    temp_fractions <- lsei::lsei(a = signatures_matrix, 
                                 b = in_mutation_catalogue_df[,i],
                                 e=G, f=H)
    temp_exposures_vector <- round(temp_fractions,digits = 6)
    names(temp_exposures_vector) <- names(in_signatures_df)
    rel_exposures_vector <- temp_exposures_vector/sum(temp_exposures_vector)
    deselect_ind <- which(rel_exposures_vector<in_per_sample_cutoff)
    temp_exposures_vector[deselect_ind] <- 0
    out_exposures_df[seq(1,dim(signatures_matrix)[2],1),i] <- 
      temp_exposures_vector
    rm(temp_fractions)
  }
  colnames(out_exposures_df) <- colnames(in_mutation_catalogue_df)
  rownames(out_exposures_df) <- colnames(in_signatures_df)
  return(out_exposures_df) 
}


#' @importFrom lsei lsei
#' 
LCD_cutoff <- function(in_mutation_catalogue_df,in_signatures_df,
                       in_cutoff=0.01,in_filename=NULL,
                       in_method="abs",in_convention="weak",
                       in_per_sample_cutoff=0) {
  # first run analysis without cutoff
  if(in_convention=="strict"){
    all_exposures_df <- LCD(in_mutation_catalogue_df,in_signatures_df)
  } else{
    all_exposures_df <- LCD(in_mutation_catalogue_df,
                            in_signatures_df,
                            in_per_sample_cutoff=in_per_sample_cutoff)
  }
  # now apply cutoff criteria to choose the main signatures
  if(in_method=="rel"){
    rel_all_exposures_df <- normalize_df_per_dim(all_exposures_df,2)
    average_rel_exposure_vector <- average_over_present(rel_all_exposures_df,1)
    sig_choice_ind <- which(average_rel_exposure_vector >= in_cutoff)
  } else {
    all_exposures_sum_df <- data.frame(sum=apply(all_exposures_df,1,sum))
    all_exposures_sum_df$sum_norm <- 
      all_exposures_sum_df$sum/sum(all_exposures_sum_df$sum)
    sig_choice_ind <- which(all_exposures_sum_df$sum_norm >= in_cutoff)
    if (!is.null(in_filename)) {
      break_vector <- seq(0,0.5,0.01)
      png(in_filename,width=400,height=400)
      hist(all_exposures_sum_df$sum_norm,breaks=break_vector,
           xlab="exposures per sig",main="Sum over whole cohort")    
      dev.off()
    }
  }
  choice_signatures_df <- in_signatures_df[,sig_choice_ind,drop=FALSE]
  # now rerun the decomposition with only the chosen signatures
  if(in_convention=="strict"){
    #out_exposures_df <- LCD_strict(in_mutation_catalogue_df,
    #                               choice_signatures_df)    
    out_exposures_df <- LCD(in_mutation_catalogue_df,choice_signatures_df)    
  } else{
    out_exposures_df <- LCD(in_mutation_catalogue_df,
                            choice_signatures_df,
                            in_per_sample_cutoff=in_per_sample_cutoff)
  }
  out_exposures_sum_df <- data.frame(sum=apply(out_exposures_df,1,sum))
  out_exposures_sum_df$sum_norm <- 
    out_exposures_sum_df$sum/sum(out_exposures_sum_df$sum)
  exposure_order <- order(out_exposures_sum_df$sum,decreasing=TRUE)
  # compute QC and error measures
  fit_catalogue_df <- as.data.frame(
    as.matrix(choice_signatures_df) %*% as.matrix(out_exposures_df))
  residual_catalogue_df <- in_mutation_catalogue_df - fit_catalogue_df
  rss <- sum(residual_catalogue_df^2)
  cosDist_fit_orig_per_matrix <- cosineDist(in_mutation_catalogue_df,
                                            fit_catalogue_df)
  cosDist_fit_orig_per_col <- rep(0,dim(in_mutation_catalogue_df)[2])
  for(i in dim(in_mutation_catalogue_df)[2]){
    cosDist_fit_orig_per_col[i] <- cosineDist(in_mutation_catalogue_df[,i],
                                              fit_catalogue_df[,i])
  }
  total_counts <- colSums(in_mutation_catalogue_df)
  sum_ind <- rev(order(total_counts))
  return(list(exposures=out_exposures_df,
              signatures=choice_signatures_df,
              choice=sig_choice_ind,
              order=exposure_order,
              residual_catalogue=residual_catalogue_df,
              rss=rss,
              cosDist_fit_orig_per_matrix=cosDist_fit_orig_per_matrix,
              cosDist_fit_orig_per_col=cosDist_fit_orig_per_col,
              sum_ind=sum_ind))
}


#'LCD with a signature-specific cutoff on exposures
#'
#'\code{LCD_cutoff} performs a mutational signatures decomposition by Linear
#'Combination Decomposition (LCD) of a given mutational catalogue \code{V} with
#'known signatures \code{W} by solving the minimization problem \eqn{min(||W*H -
#'V||)} with additional constraints of non-negativity on H where W and V are
#'known, but excludes signatures with an overall contribution less than a given
#'signature-specific cutoff (and thereby accounting for a background model) over
#'the whole cohort.
#'
#'@param in_mutation_catalogue_df A numeric data frame \code{V} with \code{n}
#'  rows and \code{m} columns, \code{n} being the number of features and
#'  \code{m} being the number of samples
#'@param in_signatures_df A numeric data frame \code{W} with \code{n} rows and
#'  \code{l} columns, \code{n} being the number of features and \code{l} being
#'  the number of signatures
#'@param in_cutoff_vector A numeric vector of values less than 1. Signatures
#'  from within \code{W} with an overall exposure less than the respective value
#'  in \code{in_cutoff_vector} will be discarded.
#'@param in_filename A path to generate a histogram of the signature exposures
#'  if non-NULL
#'@param in_method Indicate to which data the cutoff shall be applied: absolute
#'  exposures, relative exposures
#'@param in_per_sample_cutoff A numeric value less than 1. Signatures from
#'  within \code{W} with an exposure per sample less than \code{in_cutoff} will
#'  be discarded.
#'@param in_rescale Boolean, if TRUE (default) the exposures are rescaled such
#'  that colSums over exposures match colSums over mutational catalogue
#'@param in_sig_ind_df Data frame of type signature_indices_df, i.e. indicating
#'  name, function and meta-information of the signatures. Default is NULL.
#'@param in_cat_list List of categories for aggregation. Have to be among the
#'  column names of \code{in_sig_ind_df}. Default is NULL.
#'
#'@return A list with entries: \itemize{ \item \code{exposures}: The exposures
#'  \code{H}, a numeric data frame with \code{l} rows and \code{m} columns,
#'  \code{l} being the number of signatures and \code{m} being the number of
#'  samples \item \code{norm_exposures}: The normalized exposures \code{H}, a
#'  numeric data frame with \code{l} rows and \code{m} columns, \code{l} being
#'  the number of signatures and \code{m} being the number of samples \item
#'  \code{signatures}: The reduced signatures that have exposures bigger than
#'  \code{in_cutoff} \item \code{choice}: Index vector of the reduced signatures
#'  in the input signatures \item \code{order}: Order vector of the signatures
#'  by exposure \item \code{residual_catalogue}: Numerical data frame (matrix)
#'  of the difference between fit (product of signatures and exposures) and
#'  input mutational catalogue \item \code{rss}: Residual sum of squares (i.e.
#'  sum of squares of the residual catalogue) \item
#'  \code{cosDist_fit_orig_per_matrix}: Cosine distance between the fit (product
#'  of signatures and exposures) and input mutational catalogue computed after
#'  putting the matrix into vector format (i.e. one scaler product for the whole
#'  matrix) \item \code{cosDist_fit_orig_per_col}: Cosine distance between the
#'  fit (product of signatures and exposures) and input mutational catalogue
#'  computed per column (i.e. per sample, i.e. as many scaler products as there
#'  are samples in the cohort) \item \code{sum_ind}: Decreasing order of
#'  mutational loads based on the input mutational catalogue \item
#'  \code{out_sig_ind}: Data frame of the type \code{signature_indices_df}, i.e.
#'  indicating name, function and meta-information of the signatures. Default is
#'  NULL, non-NULL only if \code{in_sig_ind_df} is non-NULL. \item
#'  \code{aggregate_exposures_list}: List of exposure data frames aggregated
#'  over different categories. Default is NULL, non-NULL only if
#'  \code{in_sig_ind_df} and \code{in_cat_list} are non-NULL and if the
#'  categories specified in \code{in_cat_list} are among the column names of
#'  \code{in_sig_ind_df}. }
#'
#'@seealso \code{\link[YAPSA]{LCD}}
#'@seealso \code{\link[YAPSA]{aggregate_exposures_by_category}}
#'@seealso \code{\link[lsei]{lsei}}
#'
#' @examples
#'  NULL
#'
#'@importFrom lsei lsei
#'@export
#'
LCD_complex_cutoff <- function(in_mutation_catalogue_df,
                               in_signatures_df,
                               in_cutoff_vector=NULL,
                               in_filename=NULL,
                               in_method="abs",
                               in_per_sample_cutoff=0,
                               in_rescale=TRUE,
                               in_sig_ind_df=NULL,
                               in_cat_list=NULL) {
  # first run analysis without cutoff
  all_exposures_df <- LCD(in_mutation_catalogue_df,
                          in_signatures_df,
                          in_per_sample_cutoff=in_per_sample_cutoff)
  # now apply cutoff criteria to choose the main signatures
  if(in_method=="rel"){
    rel_all_exposures_df <- normalize_df_per_dim(all_exposures_df,2)
    average_rel_exposure_vector <- average_over_present(rel_all_exposures_df,1)
    sig_choice_ind <- which(average_rel_exposure_vector >= in_cutoff_vector &
                              average_rel_exposure_vector > 0)
    if (length(sig_choice_ind) == 0) { 
      #case that can occure in the new set of sigantures 
      return()
    }
  } else {
    all_exposures_sum_df <- data.frame(sum=apply(all_exposures_df,1,sum))
    all_exposures_sum_df$sum_norm <- 
      all_exposures_sum_df$sum/sum(all_exposures_sum_df$sum)
    sig_choice_ind <- which(all_exposures_sum_df$sum_norm >= in_cutoff_vector &
                              all_exposures_sum_df$sum_norm > 0)
    if (length(sig_choice_ind) == 0) {
      #case that can occure in the new set of sigantures 
      return()
    }
    if (!is.null(in_filename)) {
      break_vector <- seq(0,0.5,0.01)
      png(in_filename,width=400,height=400)
      hist(all_exposures_sum_df$sum_norm,breaks=break_vector,
           xlab="exposures per sig",main="Sum over whole cohort")    
      dev.off()
    }
  }
  choice_signatures_df <- in_signatures_df[,sig_choice_ind,drop=FALSE]
  # now rerun the decomposition with only the chosen signatures
  out_exposures_df <- LCD(in_mutation_catalogue_df,
                          choice_signatures_df,
                          in_per_sample_cutoff=in_per_sample_cutoff)
  out_norm_exposures_df <- normalize_df_per_dim(out_exposures_df,2)
  total_counts <- colSums(in_mutation_catalogue_df)
  sum_ind <- rev(order(total_counts))
  if(in_rescale){
    out_exposures_df <- as.data.frame(t(t(out_norm_exposures_df)*total_counts))
  }
  out_exposures_sum_df <- data.frame(sum=apply(out_exposures_df,1,sum))
  out_exposures_sum_df$sum_norm <-
    out_exposures_sum_df$sum/sum(out_exposures_sum_df$sum)
  exposure_order <- order(out_exposures_sum_df$sum,decreasing=TRUE)
  # compute QC and error measures
  fit_catalogue_df <- as.data.frame(
    as.matrix(choice_signatures_df) %*% as.matrix(out_exposures_df))
  residual_catalogue_df <- in_mutation_catalogue_df - fit_catalogue_df
  rss <- sum(residual_catalogue_df^2)
  cosDist_fit_orig_per_matrix <- cosineDist(in_mutation_catalogue_df,
                                            fit_catalogue_df)
  cosDist_fit_orig_per_col <- rep(0,dim(in_mutation_catalogue_df)[2])
  for(i in dim(in_mutation_catalogue_df)[2]){
    cosDist_fit_orig_per_col[i] <- cosineDist(in_mutation_catalogue_df[,i],
                                              fit_catalogue_df[,i])
  }
  out_sig_ind_df <- NULL
  aggregate_exposures_list <- NULL
  if(!is.null(in_sig_ind_df)){
    out_sig_ind_df <- in_sig_ind_df[sig_choice_ind,]
    if(!is.null(in_cat_list)){
      aggregate_exposures_list <- lapply(
        in_cat_list,FUN=function(current_category){
          aggregate_exposures_by_category(
            out_exposures_df,out_sig_ind_df,current_category)
        })   
      names(aggregate_exposures_list) <- in_cat_list
    }
  }
  return(list(exposures=out_exposures_df,
              norm_exposures=out_norm_exposures_df,
              signatures=choice_signatures_df,
              choice=sig_choice_ind,
              order=exposure_order,
              residual_catalogue=residual_catalogue_df,
              rss=rss,
              cosDist_fit_orig_per_matrix=cosDist_fit_orig_per_matrix,
              cosDist_fit_orig_per_col=cosDist_fit_orig_per_col,
              sum_ind=sum_ind,
              out_sig_ind_df=out_sig_ind_df,
              aggregate_exposures_list=aggregate_exposures_list))
}

#' LCD, signature-specific cutoff on exposures, per PID
#'
#'\code{\link[YAPSA]{LCD_complex_cutoff_perPID}} is a wrapper for
#'\code{\link[YAPSA]{LCD_complex_cutoff}} and runs individually for every PID.
#'
#' @export
#' @rdname LCD_complex_cutoff
#' 
LCD_complex_cutoff_perPID <- function(in_mutation_catalogue_df,
                                      in_signatures_df,
                                      in_cutoff_vector=NULL,
                                      in_filename=NULL,
                                      in_method="abs",
                                      in_rescale=TRUE,
                                      in_sig_ind_df=NULL,
                                      in_cat_list=NULL){
  complex_COSMIC_list_list <- lapply(
    seq_along(in_mutation_catalogue_df), FUN=function(current_col){
      current_mut_cat <- in_mutation_catalogue_df[,current_col,drop=FALSE]
      complex_COSMIC_list <- LCD_complex_cutoff(
        current_mut_cat,
        in_signatures_df,
        in_cutoff_vector=in_cutoff_vector,
        in_filename=NULL,
        in_method=in_method,
        in_rescale=in_rescale)
      return(complex_COSMIC_list)
    })
  exposures_list <- lapply(complex_COSMIC_list_list,FUN=function(x) {
    return(x$exposures)})
  exposures_df <- merge_exposures(exposures_list,
                                  in_signatures_df)
  names_del <- setdiff(colnames(in_mutation_catalogue_df), 
                       colnames(exposures_df)) 
                       #new makes sure that names are the same 
  if(length(names_del)>0){
    temp_df<- data.frame(replicate(length(names_del),
                                   rep(0,dim(exposures_df)[1])))
    colnames(temp_df)<-names_del
    exposures_full_df <- cbind(exposures_df, temp_df)
  }else{
    exposures_full_df<-exposures_df
  }
  norm_exposures_df <- normalize_df_per_dim(exposures_full_df,2)
  
  sig_choice_ind <- match(rownames(exposures_full_df),names(in_signatures_df))
  choice_signatures_df <- in_signatures_df[,sig_choice_ind,drop=FALSE]
  exposure_order <- order(rowSums(exposures_full_df),decreasing=TRUE)
  fit_catalogue_df <- as.data.frame(
    as.matrix(choice_signatures_df) %*% as.matrix(exposures_full_df))
  residual_catalogue_df <- in_mutation_catalogue_df - fit_catalogue_df
  rss <- sum(residual_catalogue_df^2)
  cosDist_fit_orig_per_matrix <- cosineDist(in_mutation_catalogue_df,
                                            fit_catalogue_df)
  cosDist_fit_orig_per_col <- rep(0,dim(in_mutation_catalogue_df)[2])
  for(i in dim(in_mutation_catalogue_df)[2]){
    cosDist_fit_orig_per_col[i] <- cosineDist(in_mutation_catalogue_df[,i],
                                              fit_catalogue_df[,i])
  }
  total_counts <- colSums(in_mutation_catalogue_df)
  sum_ind <- rev(order(total_counts))
  out_sig_ind_df <- NULL
  aggregate_exposures_list <- NULL
  if(!is.null(in_sig_ind_df)){
    out_sig_ind_df <- in_sig_ind_df[sig_choice_ind,]
    if(!is.null(in_cat_list)){
      aggregate_exposures_list <- lapply(
        in_cat_list,FUN=function(current_category){
          aggregate_exposures_by_category(
            exposures_full_df,out_sig_ind_df,current_category)
        })   
      names(aggregate_exposures_list) <- in_cat_list
    }
  }
  return(list(exposures=exposures_full_df,
              norm_exposures=norm_exposures_df,
              signatures=choice_signatures_df,
              choice=sig_choice_ind,
              order=exposure_order,
              residual_catalogue=residual_catalogue_df,
              rss=rss,
              cosDist_fit_orig_per_matrix=cosDist_fit_orig_per_matrix,
              cosDist_fit_orig_per_col=cosDist_fit_orig_per_col,
              sum_ind=sum_ind,
              out_sig_ind_df=out_sig_ind_df,
              aggregate_exposures_list=aggregate_exposures_list))
}


#'LCD, consensus between cohort-wide AND per-PID extraction
#'
#'\code{\link[YAPSA]{LCD_complex_cutoff_consensus}} calls
#'\code{\link[YAPSA]{LCD_complex_cutoff_combined}} AND
#'\code{\link[YAPSA]{LCD_complex_cutoff_perPID}} and makes a consensus
#'signature call set.
#'
#'@param in_cohort_LCDlist Optional, if not provided, the cohort-wide exposures
#'  are recalculated by calling \code{\link[YAPSA]{LCD_complex_cutoff}}
#'@param in_perPID_LCDlist Optional, if not provided, the per sample exposures
#'  are recalculated by calling \code{\link[YAPSA]{LCD_complex_cutoff_perPID}}
#'@param addSigs_cohort_cutoff Numeric value for a cutoff: signatures which are
#'  detected in a fraction of the samples of the cohort greater than this cutoff
#'  are kept for the consensus set of signatures
#'@param addSigs_perPID_cutoff Numeric value for a cutoff: signatures which are
#'  detected in one sample with exposure greater than this cutoff are kept for
#'  the consensus set of signatures
#'@param addSigs_relAbs_cutoff Numeric value for a cutoff: signatures which are
#'  detected with at least this fraction of all variants cohort wide are kept
#'  for the consensus set of signatures
#'@param keep.unassigned Boolean, if TRUE the exposures from the signatures
#'  which don't fulfill the criteria to be kept will be added and stored in the
#'  exposures as "unassigned", otherwise the exposures are rescaled.
#'@param keep.all.cohort.sigs If TRUE (default), all signatures extracted cohort
#'  wide are kept, if FALSE, the function reevaluates whether the signatures
#'  extracted cohort wide still fulfill their criteria (i.e. exposures > cutoff)
#'  after perPID extraction.
#'@param in_verbose Verbose if \code{in_verbose=1}
#'
#'@export
#'@rdname LCD_complex_cutoff
#'  
LCD_complex_cutoff_consensus <- function(in_mutation_catalogue_df = NULL,
                                         in_signatures_df = NULL,
                                         in_cutoff_vector = NULL,
                                         in_filename = NULL,
                                         in_method = "abs",
                                         in_rescale = TRUE,
                                         in_sig_ind_df = NULL,
                                         in_cat_list = NULL,
                                         in_cohort_LCDlist = NULL,
                                         in_perPID_LCDlist = NULL,
                                         addSigs_cohort_cutoff = 0.25,
                                         addSigs_perPID_cutoff = 0.25,
                                         addSigs_relAbs_cutoff = 0.01,
                                         keep.unassigned = FALSE,
                                         keep.all.cohort.sigs = TRUE,
                                         in_verbose = FALSE){
  ## run cohort-wide analysis
  if(is.null(in_cohort_LCDlist)){
    if(!is.null(in_mutation_catalogue_df) & !is.null(in_signatures_df)){
      in_cohort_LCDlist <- LCD_complex_cutoff(
        in_mutation_catalogue_df,
        in_signatures_df,
        in_cutoff_vector = in_cutoff_vector,
        in_filename = in_filename,
        in_method = in_method,
        in_rescale = in_rescale,
        in_sig_ind_df = in_sig_ind_df,
        in_cat_list = in_cat_list) 
    } else {
      if(in_verbose) cat("YAPSA:::LCD_complex_cutoff_combined::error:",
                         "neither in_cohort_LCDlist nor sufficient input ",
                         "provided. Abort. \n")
      return(NULL)
    }
  }
  ## run per PID analysis
  if(is.null(in_perPID_LCDlist)){
    if(!is.null(in_mutation_catalogue_df) & !is.null(in_signatures_df)){
      in_perPID_LCDlist <- LCD_complex_cutoff_perPID(
        in_mutation_catalogue_df,
        in_signatures_df,
        in_cutoff_vector = in_cutoff_vector,
        in_filename = in_filename,
        in_method = in_method,
        in_rescale = in_rescale,
        in_sig_ind_df = in_sig_ind_df,
        in_cat_list = in_cat_list) 
    } else {
      if(in_verbose) cat("YAPSA:::LCD_complex_cutoff_combined::error:",
                         "neither in_perPID_LCDlist nor sufficient input ",
                         "provided. Abort. \n")
      return(NULL)
    }
  }
  ## initialize output
  out_LCDlist <- in_perPID_LCDlist
  ## determine signatures to be removed
  if(!keep.all.cohort.sigs){
    cohortMatch_sigs <- intersect(in_cohort_LCDlist$out_sig_ind_df$sig,
                                  in_perPID_LCDlist$out_sig_ind_df$sig)
    cohortMatch_ind <- match(cohortMatch_sigs, in_sig_ind_df$sig)
    cohortMatch_CumulExpo_abs <- 
      rowSums(in_perPID_LCDlist$exposures[cohortMatch_sigs,])
    cohortMatch_CumulExpo_rel <- 
      cohortMatch_CumulExpo_abs / sum(in_perPID_LCDlist$exposures)
    select_ind <- cohortMatch_CumulExpo_rel >= in_cutoff_vector[cohortMatch_ind]
    keep_sigs <- cohortMatch_sigs[select_ind]
    keep_ind <- sort(match(keep_sigs, names(in_signatures_df)))
  }
  ## determine signatures to be added
  absExposures_sumVector <- rowSums(in_perPID_LCDlist$exposures)
  relAbsExposures_sumVector <- 
    absExposures_sumVector/sum(absExposures_sumVector)
  # addSigs_relAbs <- 
  #   names(which(relAbsExposures_sumVector > 
  #            median(relAbsExposures_sumVector)))
  addSigs_relAbs <- 
    names(which(relAbsExposures_sumVector > addSigs_relAbs_cutoff))
  # absExposures_maxVector <- apply(in_perPID_LCDlist$exposures, 1, max)
  normExposures_maxVector <- 
    apply(in_perPID_LCDlist$norm_exposures, 1, max)
  temp_df <- in_perPID_LCDlist$exposures
  temp_df[temp_df > 0] <- 1
  affected_PIDnumVector <- apply(temp_df, 1, sum)
  rel_affected_PIDnumVector <- 
    affected_PIDnumVector/dim(in_perPID_LCDlist$exposures)[2]
  addSigs_cohort <- 
    names(which(rel_affected_PIDnumVector > addSigs_cohort_cutoff))
  addSigs_perPID <- 
    names(which(normExposures_maxVector > addSigs_perPID_cutoff))
  addSigs <- 
    Reduce(union, list(addSigs_relAbs, addSigs_cohort, addSigs_perPID))
  add_ind <- sort(match(addSigs, names(in_signatures_df)))
  # relExposures_maxVector <- absExposures_maxVector/
  #    sum(absExposures_sumVector)
  # addSigs <- names(which(relExposures_maxVector > 
  #    median(relExposures_maxVector)))
  # addSigs <- names(which(normExposures_maxVector > 
  #    median(normExposures_maxVector)))
  # add_ind <- match(addSigs, names(in_signatures_df))
  # combined_cutoff <- 0.015
  # temp_cutoff_vector <- 
  #   sapply(CosmicValid_cutoffAbs_vector[names(relExposures_sumVector)],
  #          min, combined_cutoff)
  # choice_ind <- which(relExposures_sumVector >= temp_cutoff_vector)
  if(!keep.all.cohort.sigs){
    choice_ind <- sort(unique(c(keep_ind, add_ind)))
  } else {
    choice_ind <- sort(unique(c(in_cohort_LCDlist$choice, add_ind)))
  }
  out_LCDlist$choice <- choice_ind
  choice_sigs <- names(in_signatures_df)[choice_ind]
  # sig_choice_ind <- which(AlexCosmicValid_sigInd_df$sig %in% 
  #    names(choice_ind))
  # in_sig_ind_df <- AlexCosmicValid_sigInd_df[sig_choice_ind,]
  # current_cutoff_vector <- CosmicValid_cutoffAbs_vector[sig_choice_ind]
  # current_caption <- paste0("Signatures with cohort-wide exposures > ",
  #                           "combined_cutoff")
  # run_chunk(chunk_LCD_complex_cutoff_perPID)
  # out_LCDlist <- current_LCDlist
  ##  adapt output
  choice_exposures_df <- in_perPID_LCDlist$exposures[choice_sigs,]
  if(keep.unassigned){
    unassigned_vector <- colSums(in_perPID_LCDlist$exposures) - 
      colSums(choice_exposures_df)
    choice_exposures_df <- rbind(choice_exposures_df,unassigned_vector)
    rownames(choice_exposures_df)[dim(choice_exposures_df)[1]] <- "unassigned"
    out_LCDlist$exposures <- choice_exposures_df
    out_LCDlist$norm_exposures <- 
      normalize_df_per_dim(out_LCDlist$exposures, 2)
    out_LCDlist$out_sig_ind_df <- 
      in_sig_ind_df[c(choice_ind, dim(in_sig_ind_df)[1]),]
  } else {
    out_LCDlist$norm_exposures <- normalize_df_per_dim(choice_exposures_df, 2)
    total_counts <- colSums(in_mutation_catalogue_df)
    out_LCDlist$exposures <- 
      as.data.frame(t(t(out_LCDlist$norm_exposures)*total_counts))
    out_LCDlist$out_sig_ind_df <- in_sig_ind_df[choice_ind,]
  }
  out_LCDlist$signatures <- 
    in_signatures_df[,choice_sigs]
  return(out_LCDlist)
}


#' LCD, wrapper for cohort-wide, per-PID AND consensus extraction
#'
#'\code{\link[YAPSA]{LCD_complex_cutoff_combined}} is a wrapper for
#'\code{\link[YAPSA]{LCD_complex_cutoff}},
#'\code{\link[YAPSA]{LCD_complex_cutoff_perPID}} AND
#'\code{\link[YAPSA]{LCD_complex_cutoff_consensus}}.
#'
#' @export
#' @rdname LCD_complex_cutoff
#' 
LCD_complex_cutoff_combined <- function(in_mutation_catalogue_df = NULL,
                                        in_signatures_df = NULL,
                                        in_cutoff_vector = NULL,
                                        in_filename = NULL,
                                        in_method = "abs",
                                        in_rescale = TRUE,
                                        in_sig_ind_df = NULL,
                                        in_cat_list = NULL,
                                        addSigs_cohort_cutoff = 0.25,
                                        addSigs_perPID_cutoff = 0.25,
                                        addSigs_relAbs_cutoff = 0.01,
                                        keep.all.cohort.sigs = TRUE,
                                        in_verbose = FALSE){
  cohort_LCDlist <- LCD_complex_cutoff(
    in_mutation_catalogue_df,
    in_signatures_df,
    in_cutoff_vector = in_cutoff_vector,
    in_filename = in_filename,
    in_method = in_method,
    in_rescale = in_rescale,
    in_sig_ind_df = in_sig_ind_df,
    in_cat_list = in_cat_list)
  perPID_LCDlist <- LCD_complex_cutoff_perPID(
    in_mutation_catalogue_df,
    in_signatures_df,
    in_cutoff_vector = in_cutoff_vector,
    in_filename = in_filename,
    in_method = in_method,
    in_rescale = in_rescale,
    in_sig_ind_df = in_sig_ind_df,
    in_cat_list = in_cat_list) 
  consensus_LCDlist <- LCD_complex_cutoff_consensus(
    in_mutation_catalogue_df = in_mutation_catalogue_df,
    in_signatures_df = in_signatures_df,
    in_sig_ind_df = in_sig_ind_df,
    in_cohort_LCDlist = cohort_LCDlist,
    in_perPID_LCDlist = perPID_LCDlist,
    addSigs_cohort_cutoff = addSigs_cohort_cutoff,
    addSigs_perPID_cutoff = addSigs_perPID_cutoff,
    addSigs_relAbs_cutoff = addSigs_relAbs_cutoff,
    keep.unassigned = TRUE,
    keep.all.cohort.sigs = keep.all.cohort.sigs)
  consensusRescale_LCDlist <- LCD_complex_cutoff_consensus(
    in_mutation_catalogue_df = in_mutation_catalogue_df,
    in_signatures_df = in_signatures_df,
    in_sig_ind_df = in_sig_ind_df,
    in_cohort_LCDlist = cohort_LCDlist,
    in_perPID_LCDlist = perPID_LCDlist,
    addSigs_cohort_cutoff = addSigs_cohort_cutoff,
    addSigs_perPID_cutoff = addSigs_perPID_cutoff,
    addSigs_relAbs_cutoff = addSigs_relAbs_cutoff,
    keep.unassigned = FALSE,
    keep.all.cohort.sigs = keep.all.cohort.sigs)
  return(list(cohort = cohort_LCDlist,
              perPID = perPID_LCDlist,
              consensus = consensus_LCDlist,
              consensusRescale = consensusRescale_LCDlist))
}


res <- function(x,b,in_matrix){
  b - in_matrix %*% x
}


norm_res <- function(x,b,in_matrix){
  norm(as.matrix(b - in_matrix %*% x),"F")
}



#' CD stratification analysis
#'
#' @param in_mutation_sub_catalogue_list A list of \code{s} stratified
#'   mutational catalogues \code{Vi} \(numeric data frames\) with \code{n} rows
#'   and \code{m} columns each, \code{n} being the number of features and
#'   \code{m} being the number of samples. This list is naturally provided in
#'   \code{\link{run_SMC}}.
#' @param in_signatures_df A numeric data frame \code{W} with \code{n} rows and
#'   \code{l} columns, \code{n} being the number of features and \code{l} being
#'   the number of signatures
#' @param in_F_df Default NULL
#'
#' @return Returns a list with all exposures and the stratified ones
#'
#' @importFrom lsei lsei
#'   
LCD_SMC <- function(in_mutation_sub_catalogue_list,
                    in_signatures_df,in_F_df=NULL){
  ## find general properties
  number_of_strata <- length(in_mutation_sub_catalogue_list)
  number_of_sigs <- dim(in_signatures_df)[2]
  number_of_PIDs <- dim(in_mutation_sub_catalogue_list[[1]])[2]
  number_of_features <- dim(in_mutation_sub_catalogue_list[[1]])[1]
  ## 1. construct composite signatures_matrix
  signatures_matrix_element <- as.matrix(in_signatures_df)
  zero_element <- matrix(
    rep(0,dim(signatures_matrix_element)[1]*dim(signatures_matrix_element)[2]),
    ncol=dim(signatures_matrix_element)[2])
  temp_matrix <- NULL
  for (i in seq_len(number_of_strata)) {
    temp_row <- NULL
    for (j in seq_len(number_of_strata)) {
      if (i==j) {
        temp_row <- cbind(temp_row,signatures_matrix_element)
      } else {
        temp_row <- cbind(temp_row,zero_element)
      }
    }
    temp_matrix <- rbind(temp_matrix,temp_row)
  }
  signatures_matrix <- temp_matrix
  ## 2. construct composite mutation catalogue
  temp_matrix <- NULL
  for (i in seq_len(number_of_strata)) {
    temp_matrix <- rbind(temp_matrix,in_mutation_sub_catalogue_list[[i]])
  }
  pasted_mutation_catalogue_df <- temp_matrix
  ## 3. account for boundary conditions
  ## 3.a) account for equality boundary condition
  if (!is.null(in_F_df)) {
    # This condition is fulfilled when an exposures file has been supplied.
    F_df <- in_F_df
  } else {
    # This is the standard case when no exposures file has been supplied and 
    # the exposures have to be computed by LCD.
    sum_df <- data.frame(matrix(rep(0,number_of_PIDs*number_of_features),
                                ncol=number_of_PIDs))
    for (i in seq_len(number_of_strata)) {
      sum_df <- sum_df + in_mutation_sub_catalogue_list[[i]]
    }
    all_mutation_catalogue_df <- sum_df
    F_df <- LCD(all_mutation_catalogue_df,in_signatures_df)    
  }
  diagonal_element <- diag(number_of_sigs)
  temp_matrix <- NULL
  for (i in seq_len(number_of_strata)) {
    temp_matrix <- cbind(temp_matrix,diagonal_element)
  }
  E <- temp_matrix
  ## 3.b) account for inequality boundary condition
  G <- diag(dim(signatures_matrix)[2])
  H <- rep(0,dim(signatures_matrix)[2])
  out_exposures_df <- data.frame()
  for (i in seq_len(ncol(pasted_mutation_catalogue_df))) {
    # temp_fractions <- limSolve::lsei(A = signatures_matrix, 
    #                                  B = pasted_mutation_catalogue_df[,i], 
    #                                  E=E, F=F_df[,i], G=G, H=H)
    temp_fractions <- lsei::lsei(a = signatures_matrix, 
                                 b = pasted_mutation_catalogue_df[,i], 
                                 c=E, d=F_df[,i], e=G, f=H)
    temp_exposures_vector <- round(temp_fractions,digits = 6)
    names(temp_exposures_vector) <- names(in_signatures_df)
    out_exposures_df[seq(1,dim(signatures_matrix)[2],1),i] <- 
      #as.vector(temp_fractions$X)
      as.vector(temp_exposures_vector)
    rm(temp_fractions)
  }
  out_list <- list()
  for (i in seq_len(number_of_strata)) {
    out_list[[i]] <- as.data.frame(
      out_exposures_df[seq((number_of_sigs*(i-1)+1),(number_of_sigs*i),1),])
    colnames(out_list[[i]]) <- colnames(in_mutation_sub_catalogue_list[[1]])
    rownames(out_list[[i]]) <- colnames(in_signatures_df)
  }
  colnames(F_df) <- colnames(in_mutation_sub_catalogue_list[[1]])
  rownames(F_df) <- colnames(in_signatures_df)
  return(list(exposures_all_df=F_df,sub_exposures_list=out_list))
}
