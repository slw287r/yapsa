## ----load_style, warning=FALSE, echo=FALSE, message=FALSE, results="hide"-----
library(BiocStyle)

## ----packages, include=FALSE--------------------------------------------------
library(YAPSA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(knitr)
opts_chunk$set(echo=TRUE)
opts_chunk$set(fig.show='asis')

## ---- load_stored_sig_data----------------------------------------------------
data(sigs)
data(cutoffs)
data("lymphomaNature2013_mutCat_df")
current_cutoff_vector <- cutoffCosmicValid_abs_df[6,]

## ----LCD with cutoffs---------------------------------------------------------
lymphoma_COSMIC_listsList <-
  LCD_complex_cutoff_combined(
      in_mutation_catalogue_df = lymphomaNature2013_mutCat_df,
      in_cutoff_vector = current_cutoff_vector, 
      in_signatures_df = AlexCosmicValid_sig_df, 
      in_sig_ind_df = AlexCosmicValid_sigInd_df)

## ----subrgroup annotation-----------------------------------------------------
data(lymphoma_PID)
colnames(lymphoma_PID_df) <- "SUBGROUP"
lymphoma_PID_df$PID <- rownames(lymphoma_PID_df)
COSMIC_subgroups_df <- 
  make_subgroups_df(lymphoma_PID_df,
                    lymphoma_COSMIC_listsList$cohort$exposures)

## ----caption_exposures, echo=FALSE--------------------------------------------
cap <- "Exposures to SNV mutational signatures"

## ----exposures_cutoffs, warning=FALSE, fig.width=8, fig.height=6, fig.cap= cap----
exposures_barplot(
  in_exposures_df = lymphoma_COSMIC_listsList$cohort$exposures,
  in_signatures_ind_df = lymphoma_COSMIC_listsList$cohort$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)               

## ----format_raw, warning=FALSE, fig.width=8, fig.height=6, fig.cap= cap-------
data("lymphoma_Nature2013_raw")
names(lymphoma_PID_df) <- gsub("SUBGROUP", "subgroup", names(lymphoma_PID_df))
names(lymphoma_Nature2013_raw_df) <- c("PID","TYPE","CHROM","START",
                                       "STOP","REF","ALT","FLAG")
lymphoma_Nature2013_df <- subset(lymphoma_Nature2013_raw_df,TYPE=="subs",
                                 select=c(CHROM,START,REF,ALT,PID))
names(lymphoma_Nature2013_df)[2] <- "POS"
lymphoma_Nature2013_df$SUBGROUP <- "unknown"
DLBCL_ind <- grep("^DLBCL.*",lymphoma_Nature2013_df$PID)
lymphoma_Nature2013_df$SUBGROUP[DLBCL_ind] <- "DLBCL_other"
MMML_ind <- grep("^41[0-9]+$",lymphoma_Nature2013_df$PID)
lymphoma_Nature2013_df <- lymphoma_Nature2013_df[MMML_ind,]
for(my_PID in rownames(lymphoma_PID_df)) {
  PID_ind <- which(as.character(lymphoma_Nature2013_df$PID)==my_PID)
  lymphoma_Nature2013_df$SUBGROUP[PID_ind] <-
    lymphoma_PID_df$subgroup[which(rownames(lymphoma_PID_df)==my_PID)]
}
lymphoma_Nature2013_df$SUBGROUP <- factor(lymphoma_Nature2013_df$SUBGROUP)
lymphoma_Nature2013_df <- translate_to_hg19(lymphoma_Nature2013_df,"CHROM")
lymphoma_Nature2013_df$change <- 
  attribute_nucleotide_exchanges(lymphoma_Nature2013_df)
lymphoma_Nature2013_df <- 
  lymphoma_Nature2013_df[order(lymphoma_Nature2013_df$PID,
                               lymphoma_Nature2013_df$CHROM,
                               lymphoma_Nature2013_df$POS),]
lymphoma_Nature2013_df <- annotate_intermut_dist_cohort(lymphoma_Nature2013_df,
                                                        in_PID.field="PID")
data("exchange_colour_vector")
lymphoma_Nature2013_df$col <- 
  exchange_colour_vector[lymphoma_Nature2013_df$change]

## ----stratify_mut_density-----------------------------------------------------
lymphoma_Nature2013_df$density_cat <- cut(lymphoma_Nature2013_df$dist,
                                          c(0,1001,100001,Inf),
                                          right=FALSE,
                                          labels=c("high","intermediate",
                                                   "background"))

## ----table_strata_mut_density-------------------------------------------------
temp_df <- data.frame(table(lymphoma_Nature2013_df$density_cat))
names(temp_df) <- c("Stratum","Cohort-wide counts")
kable(temp_df, caption=paste0("Strata for the SMC of mutation density"))

## ----caption_SMC, echo=FALSE--------------------------------------------------
cap="SMC (Stratification of the Mutational Catalogue)
        based on mutation density."

## ----SMC, results="hide",warning=FALSE,fig.width=8,fig.height=7,fig.cap=cap----
strata_order_ind <- c(1,3,2)
mut_density_list <- run_SMC(lymphoma_Nature2013_df,
                            lymphoma_COSMIC_listsList$cohort$signatures,
                            lymphoma_COSMIC_listsList$cohort$out_sig_ind_df,
                            COSMIC_subgroups_df,
                            column_name="density_cat",
                            refGenome=BSgenome.Hsapiens.UCSC.hg19,
                            cohort_method_flag="norm_PIDs",
                            in_strata_order_ind=strata_order_ind)

## ----caption_SMC_dodged, echo=FALSE-------------------------------------------
cap = "SMC results displayed as dodged plot."

## ----SMC_dodged, results="hide",warning=FALSE,fig.width=8,fig.height=7,fig.cap=cap----
dodged_df <- do.call(rbind, mut_density_list$cohort)
names(dodged_df) <- gsub("variable","stratum", names(dodged_df))
names(dodged_df) <- gsub("sig","signature", names(dodged_df))
dodged_df$stratum <- gsub("_rel", "", as.character(dodged_df$stratum))
dodged_df <- dodged_df[which(dodged_df$stratum != "all"),]
dodged_df$stratum <-
  factor(as.character(dodged_df$stratum),
         levels = sort(unique(dodged_df$stratum))[rev(strata_order_ind)])
dodged_plot <- ggplot() +
  geom_bar(data = dodged_df,
           aes_string(x = "signature", y = "exposure",
                      group = "stratum", fill = "stratum"),
           stat = "identity", position = "dodge", size = 1.5) +
  geom_errorbar(data = dodged_df,
                aes_string(x = "signature", ymin = "exposure_min",
                           ymax = "exposure_max",
                           group = "stratum"),
                position = position_dodge(width = 0.9), width = 0.3) +
  labs(y = "relative exposures")
if(exists("current_strata_colVector")){
  dodged_plot <- dodged_plot +
    scale_fill_manual(values = current_strata_colVector[-1],
                      labels = current_labelVector[-1]
                                [current_strata_order_ind])
}
print(dodged_plot)

## ----stat_SMC_mut_density, warning=FALSE, message=FALSE-----------------------
stat_mut_density_list <- stat_test_SMC(mut_density_list,in_flag="norm")
kable(stat_mut_density_list$kruskal_df,
      caption=paste0("Results of Kruskal tests for cohort-wide exposures over",
                     " strata per signature without and with correction for ",
                     "multiple testing."))

## ----post_hoc_mut_density-----------------------------------------------------
significance_level <- 0.05
for(i in seq_len(dim(stat_mut_density_list$kruskal_df)[1])){
  if(stat_mut_density_list$kruskal_df$Kruskal_p_val_BH[i]<significance_level){
    print(paste0("Signature: ",rownames(stat_mut_density_list$kruskal_df)[i]))
    print(stat_mut_density_list$kruskal_posthoc_list[[i]])
  }
}

