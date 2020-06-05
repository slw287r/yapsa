# Copyright © 2014-2019  The YAPSA package contributors
# This file is part of the YAPSA package. The YAPSA package is licenced under
# GPL-3

#'INDEL function V1 - not compartible with AlexandrovSignatures
#' @param in_df Input data frame containing the variances in a vcf-like format 
#' @param in_ALT.field Column number for alternitve field
#' @param in_REF.field Coloumn number for reference field  
#' @param in_breaks Handed over to function cut
#' @param in_labels Handed over to function cut
#'
#' @return classVector, a factor vector of indel sizes
#' @export
#'
#' @examples
#'  NULL
#' 
classify_indels <- function(in_df,
                            in_ALT.field = "ALT",
                            in_REF.field = "REF",
                            in_breaks = c(-Inf, -10, -3, 0, 2, 9, Inf),
                            in_labels = c("del3", "del2", "del1", 
                                          "in1", "in2", "in3")){
  deltaVector <- nchar(as.character(in_df[, in_ALT.field])) - 
    nchar(as.character(in_df[, in_REF.field]))
  classVector <- cut(deltaVector, breaks = in_breaks, labels = in_labels)
  return(classVector)
}
