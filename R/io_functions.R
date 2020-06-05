# Copyright © 2014-2019  The YAPSA package contributors
# This file is part of the YAPSA package. The YAPSA package is licenced under
# GPL-3

#' Read a single vcf-like file into a single data frame
#'
#' Note: this function uses \code{\link[utils]{read.csv}} to read vcf-like files
#' into data frames for single samples. As it uses
#' \code{\link[utils]{read.csv}}, the default value for \code{comment.char} is
#' "" and not "#" as it would have been for \code{\link[utils]{read.table}}.
#'
#' @param current_ind Index of the file to read from the list provided below.
#' @param in_list List of paths to vcf-like file to be read. The list may be
#'   named.
#' @param header Boolean whether a header information should be read (as in
#'   \code{\link[utils]{read.table}})
#' @param in_header Vector of column names to be substituted if non-NULL.
#' @param variant_type Default is "SNV" and provides additional plausibility and
#'   checks, omitted if other string
#' @param delete.char Character to be deleted, e.g. in order to discriminate
#'   between comment lines and header lines, if non-NULL
#' @param ... Parameters passed on to \code{\link[utils]{read.table}}
#'
#' @examples
#'  NULL
#'
#' @return A vcf-like data frame
#' @export
#' 
read_entry <- function(current_ind,
                       in_list,
                       header = TRUE,
                       in_header = NULL,
                       variant_type = "SNV",
                       delete.char = NULL,
                       ...){
  current_PID <- names(in_list)[current_ind]
  current_entry <- in_list[[current_ind]]
  vcf_like_df <- tryCatch({
    if(is.null(delete.char)){
      temp_df <- read.csv(current_entry, header = header, ...)
    } else {
      temp_df <- read.csv(pipe(paste0('zcat ', current_entry, " | sed s/^", 
                                        delete.char, "//")), 
                            header = header, ...)
    }
    if(!header & !is.null(in_header)) names(temp_df) <- in_header
    ## remove false counts, i.e. keep only SNVs
    if(variant_type == "SNV"){
      temp_true_ind <- which((temp_df$REF %in% c("A","C","G","T")) &
                               (temp_df$ALT %in% c("A","C","G","T")))
      temp_df <- temp_df[temp_true_ind,]
    }
    ## attribute PID
    temp_df$PID <- current_PID
    temp_df
  },
  error=function(cond){
    message(paste0("read_entry::error. Original error message:"))
    message(paste0(cond))
    message(paste0("Return NULL.\n"))
    return(NULL)
  })
  return(vcf_like_df)
}


#' Read a list of vcf-like files into a list of data frames
#'
#' @param in_parallel If multicore functionality is provided on a compute
#'   cluster, this option may be set to TRUE in order to enhance speed.
#'
#' @examples
#'  NULL
#'
#' @return A list with entries: \itemize{ \item \code{vcf_like_df_list}: List of
#'   the read data frames \item \code{readVcf_time}: Object of class
#'   \code{proc_time}, which stores the time needed for reading in the data }
#'
#' @importFrom doParallel registerDoParallel
#' @export
#'
#' @rdname read_entry
#'   
read_list <- function(in_list,
                      in_parallel = FALSE,
                      header = TRUE,
                      in_header = NULL,
                      ...){
  seq_list <- seq_along(in_list)
  names(seq_list) <- names(in_list)
  buildCatalogues_time <- 0
  if(in_parallel){
    #library(parallel)
    detectCores()
    cl <- makeCluster(detectCores() - 1)
    #library(parallel)
    #library(doParallel)
    registerDoParallel(cl, cores = detectCores() - 1)
    start_time <- proc.time()
    mut_cat_df_list <- mclapply(seq_along(in_list),function(current_ind)
      read_entry(current_ind, in_list=in_list, header = header,
                 in_header = in_header, ...))
    buildCatalogues_time <- proc.time() - start_time
    stopCluster(cl)
  } else {
    start_time <- proc.time()
    vcf_like_df_list <- lapply(seq_list,function(current_ind)
      read_entry(current_ind, in_list=in_list, header = header,
                 in_header = in_header, ...))
    readVcf_time <- proc.time() - start_time
  }
  vcf_like_df_list <- vcf_like_df_list[
    which(!unlist(lapply(vcf_like_df_list,is.null)))]
  return(list(vcf_like_df_list=vcf_like_df_list,
              readVcf_time=readVcf_time))
}
