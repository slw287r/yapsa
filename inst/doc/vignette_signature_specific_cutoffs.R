## ----load_style, warning=FALSE, echo=FALSE, message=FALSE, results="hide"-----
library(BiocStyle)

## ----packages, include=FALSE--------------------------------------------------
library(YAPSA)
library(knitr)
library(gridExtra)
library(dplyr)
opts_chunk$set(echo=TRUE)
opts_chunk$set(fig.show='asis')

## ----load signatures----------------------------------------------------------
data(sigs)
data(sigs_pcawg)

## ----load cutoffs-------------------------------------------------------------
data(cutoffs)
data(cutoffs_pcawg)

## ----cutoff example-----------------------------------------------------------
data(cutoffs)
current_cutoff_vector <- cutoffCosmicValid_abs_df[6, ]
current_cutoff_vector

## ----load data----------------------------------------------------------------
data(lymphomaNature2013_mutCat_df)

## ----set variables------------------------------------------------------------
current_sig_df <- AlexCosmicValid_sig_df
current_sigInd_df <- AlexCosmicValid_sigInd_df

## ----cutoff vector------------------------------------------------------------
current_cutoff_vector <- rep(0, dim(AlexCosmicValid_sig_df)[2])

## ----lymphoma_cohort_LCD_results----------------------------------------------
lymphoma_COSMIC_zero_listsList <-
  LCD_complex_cutoff_combined(
    in_mutation_catalogue_df = lymphomaNature2013_mutCat_df,
    in_cutoff_vector = current_cutoff_vector, 
    in_signatures_df = current_sig_df, 
    in_sig_ind_df = current_sigInd_df)

## ----subrgroup annotation-----------------------------------------------------
data(lymphoma_PID)
colnames(lymphoma_PID_df) <- "SUBGROUP"
lymphoma_PID_df$PID <- rownames(lymphoma_PID_df)
COSMIC_subgroups_df <- 
  make_subgroups_df(lymphoma_PID_df,
                    lymphoma_COSMIC_zero_listsList$cohort$exposures)

## ----caption_barplot_2, echo=FALSE--------------------------------------------
cap <- "Absoute exposures of the COSMIC signatures in the lymphoma mutational
        catalogs, signature-specific cutoffs with a cost factor of 6 used
        for the LCD"

## ----exposures_zero, warning=FALSE, fig.width=8, fig.height=6, fig.cap= cap----
result_cohort <- lymphoma_COSMIC_zero_listsList$cohort
exposures_barplot(
  in_exposures_df = result_cohort$exposures,
  in_signatures_ind_df = result_cohort$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)        

## ----set signature-specific cutoff--------------------------------------------
current_cutoff_df <- cutoffCosmicValid_abs_df
current_cost_factor <- 6
current_cutoff_vector <- current_cutoff_df[current_cost_factor,]

## ----LCD with cutoffs---------------------------------------------------------
lymphoma_COSMIC_listsList <-
  LCD_complex_cutoff_combined(
      in_mutation_catalogue_df = lymphomaNature2013_mutCat_df,
      in_cutoff_vector = current_cutoff_vector, 
      in_signatures_df = current_sig_df, 
      in_sig_ind_df = current_sigInd_df)

## ----caption_barplot, echo=FALSE----------------------------------------------
cap <- "Absolute exposures of the COSMIC signatures in the lymphoma mutational
        catalogs, signature-specific cutoffs with a cost factor of 6 used
        for the LCD"

## ----exposures_cutoffs, warning=FALSE, fig.width=8, fig.height=6, fig.cap= cap----
result_cohort <- lymphoma_COSMIC_listsList$cohort
exposures_barplot(
  in_exposures_df = result_cohort$exposures,
  in_signatures_ind_df = result_cohort$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)               

