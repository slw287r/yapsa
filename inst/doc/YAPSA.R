## ----load_style, warning=FALSE, message=FALSE, results="hide"-----------------
library(BiocStyle)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(YAPSA)
library(knitr)
opts_chunk$set(echo=TRUE)
opts_chunk$set(fig.show='asis')

## ----load_lymphoma------------------------------------------------------------
data("lymphoma_Nature2013_raw")

## ----load_lymphoma_ftp, eval=FALSE--------------------------------------------
#  lymphoma_Nature2013_ftp_path <- paste0(
#    "ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/",
#    "somatic_mutation_data/Lymphoma B-cell/",
#    "Lymphoma B-cell_clean_somatic_mutations_",
#    "for_signature_analysis.txt")
#  lymphoma_Nature2013_raw_df <- read.csv(file=lymphoma_Nature2013_ftp_path,
#                                         header=FALSE,sep="\t")

## ----adapt_data---------------------------------------------------------------
names(lymphoma_Nature2013_raw_df) <- c("PID","TYPE","CHROM","START",
                                       "STOP","REF","ALT","FLAG")
lymphoma_Nature2013_df <- subset(lymphoma_Nature2013_raw_df,TYPE=="subs",
                                 select=c(CHROM,START,REF,ALT,PID))
names(lymphoma_Nature2013_df)[2] <- "POS"
kable(head(lymphoma_Nature2013_df))

## ----PID_info-----------------------------------------------------------------
unique(lymphoma_Nature2013_df$PID)

## ----annotate_subgroup, warning=FALSE, message=FALSE--------------------------
lymphoma_Nature2013_df$SUBGROUP <- "unknown"
DLBCL_ind <- grep("^DLBCL.*",lymphoma_Nature2013_df$PID)
lymphoma_Nature2013_df$SUBGROUP[DLBCL_ind] <- "DLBCL_other"
MMML_ind <- grep("^41[0-9]+$",lymphoma_Nature2013_df$PID)
lymphoma_Nature2013_df <- lymphoma_Nature2013_df[MMML_ind,]
data(lymphoma_PID)
for(my_PID in rownames(lymphoma_PID_df)) {
  PID_ind <- which(as.character(lymphoma_Nature2013_df$PID)==my_PID)
  lymphoma_Nature2013_df$SUBGROUP[PID_ind] <-
    lymphoma_PID_df$subgroup[which(rownames(lymphoma_PID_df)==my_PID)]
}
lymphoma_Nature2013_df$SUBGROUP <- factor(lymphoma_Nature2013_df$SUBGROUP)
unique(lymphoma_Nature2013_df$SUBGROUP)

## ----intermut_distances, warning=FALSE, message=FALSE, results="hide"---------
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

## ----make_rainfall_plot, fig.cap="Rainfall plot in a trellis structure."------
choice_PID <- "4121361"
PID_df <- subset(lymphoma_Nature2013_df,PID==choice_PID)
#trellis_rainfall_plot(PID_df,in_point_size=unit(0.5,"mm"))

## ---- load_stored_sig_data----------------------------------------------------
data(sigs)

## ----process_old_sigs---------------------------------------------------------
Alex_signatures_path <- paste0("ftp://ftp.sanger.ac.uk/pub/cancer/",
                               "AlexandrovEtAl/signatures.txt")
AlexInitialArtif_sig_df <- read.csv(Alex_signatures_path,header=TRUE,sep="\t")
kable(AlexInitialArtif_sig_df[c(1:9),c(1:4)])

## ----reformat_old_sigs--------------------------------------------------------
Alex_rownames <- paste(AlexInitialArtif_sig_df[,1],
                       AlexInitialArtif_sig_df[,2],sep=" ")
select_ind <- grep("Signature",names(AlexInitialArtif_sig_df))
AlexInitialArtif_sig_df <- AlexInitialArtif_sig_df[,select_ind]
number_of_Alex_sigs <- dim(AlexInitialArtif_sig_df)[2]
names(AlexInitialArtif_sig_df) <- gsub("Signature\\.","A",
                                       names(AlexInitialArtif_sig_df))
rownames(AlexInitialArtif_sig_df) <- Alex_rownames
kable(AlexInitialArtif_sig_df[c(1:9),c(1:6)],
      caption="Exemplary data from the initial Alexandrov signatures.")

## ----filter_validated_signatures----------------------------------------------
AlexInitialValid_sig_df <- AlexInitialArtif_sig_df[,grep("^A[0-9]+",
                                          names(AlexInitialArtif_sig_df))]
number_of_Alex_validated_sigs <- dim(AlexInitialValid_sig_df)[2]

## ----load_sigs----------------------------------------------------------------
Alex_COSMIC_signatures_path <- 
  paste0("http://cancer.sanger.ac.uk/cancergenome/",
         "assets/signatures_probabilities.txt")
AlexCosmicValid_sig_df <- read.csv(Alex_COSMIC_signatures_path,
                                      header=TRUE,sep="\t")
Alex_COSMIC_rownames <- paste(AlexCosmicValid_sig_df[,1],
                              AlexCosmicValid_sig_df[,2],sep=" ")
COSMIC_select_ind <- grep("Signature",names(AlexCosmicValid_sig_df))
AlexCosmicValid_sig_df <- AlexCosmicValid_sig_df[,COSMIC_select_ind]
number_of_Alex_COSMIC_sigs <- dim(AlexCosmicValid_sig_df)[2]
names(AlexCosmicValid_sig_df) <- gsub("Signature\\.","AC",
                                         names(AlexCosmicValid_sig_df))
rownames(AlexCosmicValid_sig_df) <- Alex_COSMIC_rownames
kable(AlexCosmicValid_sig_df[c(1:9),c(1:6)],
      caption="Exemplary data from the updated Alexandrov signatures.")

## ----reorder_features---------------------------------------------------------
COSMIC_order_ind <- match(Alex_rownames,Alex_COSMIC_rownames)
AlexCosmicValid_sig_df <- AlexCosmicValid_sig_df[COSMIC_order_ind,]
kable(AlexCosmicValid_sig_df[c(1:9),c(1:6)],
      caption=paste0("Exemplary data from the updated Alexandrov ",
                     "signatures, rows reordered."))

## ---- build_sig_ind_df--------------------------------------------------------
signature_colour_vector <- c("darkgreen","green","pink","goldenrod",
                             "lightblue","blue","orangered","yellow",
                             "orange","brown","purple","red",
                             "darkblue","magenta","maroon","yellowgreen",
                             "violet","lightgreen","sienna4","deeppink",
                             "darkorchid","seagreen","grey10","grey30",
                             "grey50","grey70","grey90")
bio_process_vector <- c("spontaneous deamination","spontaneous deamination",
                        "APOBEC","BRCA1_2","Smoking","unknown",
                        "defect DNA MMR","UV light exposure","unknown",
                        "IG hypermutation","POL E mutations","temozolomide",
                        "unknown","APOBEC","unknown","unknown","unknown",
                        "unknown","unknown","unknown","unknown","unknown",
                        "nonvalidated","nonvalidated","nonvalidated",
                        "nonvalidated","nonvalidated")
AlexInitialArtif_sigInd_df <- data.frame(sig=colnames(AlexInitialArtif_sig_df))
AlexInitialArtif_sigInd_df$index <- seq_len(dim(AlexInitialArtif_sigInd_df)[1])
AlexInitialArtif_sigInd_df$colour <- signature_colour_vector
AlexInitialArtif_sigInd_df$process <- bio_process_vector

COSMIC_signature_colour_vector <- c("green","pink","goldenrod",
                                    "lightblue","blue","orangered","yellow",
                                    "orange","brown","purple","red",
                                    "darkblue","magenta","maroon",
                                    "yellowgreen","violet","lightgreen",
                                    "sienna4","deeppink","darkorchid",
                                    "seagreen","grey","darkgrey",
                                    "black","yellow4","coral2","chocolate2",
                                    "navyblue","plum","springgreen")
COSMIC_bio_process_vector <- c("spontaneous deamination","APOBEC",
                               "defect DNA DSB repair hom. recomb.",
                               "tobacco mutatgens, benzo(a)pyrene",
                               "unknown",
                               "defect DNA MMR, found in MSI tumors",
                               "UV light exposure","unknown","POL eta and SHM",
                               "altered POL E",
                               "alkylating agents, temozolomide",
                               "unknown","APOBEC","unknown",
                               "defect DNA MMR","unknown","unknown",
                               "unknown","unknown",
                               "associated w. small indels at repeats",
                               "unknown","aristocholic acid","unknown",
                               "aflatoxin","unknown","defect DNA MMR",
                               "unknown","unknown","tobacco chewing","unknown")
AlexCosmicValid_sigInd_df <- data.frame(sig=colnames(AlexCosmicValid_sig_df))
AlexCosmicValid_sigInd_df$index <- seq_len(dim(AlexCosmicValid_sigInd_df)[1])
AlexCosmicValid_sigInd_df$colour <- COSMIC_signature_colour_vector
AlexCosmicValid_sigInd_df$process <- COSMIC_bio_process_vector

## ----load_libraries, warning=FALSE, message=FALSE-----------------------------
library(BSgenome.Hsapiens.UCSC.hg19)

## ----build_mutational_catalogue, results="hide"-------------------------------
word_length <- 3

lymphomaNature2013_mutCat_list <- 
  create_mutation_catalogue_from_df(
    lymphoma_Nature2013_df,
    this_seqnames.field = "CHROM", this_start.field = "POS",
    this_end.field = "POS", this_PID.field = "PID",
    this_subgroup.field = "SUBGROUP",
    this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
    this_wordLength = word_length)

## ----display_structure_mutation_catalogue_list--------------------------------
names(lymphomaNature2013_mutCat_list)
lymphomaNature2013_mutCat_df <- as.data.frame(
  lymphomaNature2013_mutCat_list$matrix)
kable(lymphomaNature2013_mutCat_df[c(1:9),c(5:10)])

## ----LCD, warning=FALSE-------------------------------------------------------
current_sig_df <- AlexCosmicValid_sig_df
current_sigInd_df <- AlexCosmicValid_sigInd_df
lymphomaNature2013_COSMICExposures_df <-
  LCD(lymphomaNature2013_mutCat_df,current_sig_df)

## ----make_subgroups_df_LCD----------------------------------------------------
COSMIC_subgroups_df <- 
  make_subgroups_df(lymphoma_Nature2013_df,
                    lymphomaNature2013_COSMICExposures_df)

## ----caption_barplot, echo=FALSE----------------------------------------------
cap <- "Absoute exposures of the COSMIC signatures in the lymphoma mutational
        catalogues, no cutoff for the LCD (Linear Combination Decomposition)"

## ----exposures_barplot_LCD, fig.width=6, fig.height=5, fig.cap=cap------------
exposures_barplot(
  in_exposures_df = lymphomaNature2013_COSMICExposures_df,
  in_subgroups_df = COSMIC_subgroups_df)

## ----caption_barplot_2, echo=FALSE--------------------------------------------
cap <- "Absoute exposures of the COSMIC signatures in the lymphoma mutational
    catalogues, no cutoff for the LCD (Linear Combination Decomposition)"

## ----exposures_barplot_LCD_sig_ind, fig.width=6, fig.height=5, fig.cap=cap----
exposures_barplot(
  in_exposures_df = lymphomaNature2013_COSMICExposures_df,
  in_signatures_ind_df = current_sigInd_df,
  in_subgroups_df = COSMIC_subgroups_df)

## ----rowSums_in_cohort--------------------------------------------------------
rowSums(lymphomaNature2013_COSMICExposures_df)

## ----LCD_complex_cutoff_cutoffZero--------------------------------------------
zero_cutoff_vector <- rep(0,dim(current_sig_df)[2])
CosmicValid_cutoffZero_LCDlist <- LCD_complex_cutoff(
  in_mutation_catalogue_df = lymphomaNature2013_mutCat_df,
  in_signatures_df = current_sig_df,
  in_cutoff_vector = zero_cutoff_vector,
  in_sig_ind_df = current_sigInd_df)

## ----apply_make_subgroups_df_cutoffZero---------------------------------------
COSMIC_subgroups_df <- 
  make_subgroups_df(lymphoma_Nature2013_df,
                    CosmicValid_cutoffZero_LCDlist$exposures)

## ----caption_barplot_3, echo=FALSE--------------------------------------------
cap="Absoute exposures of the COSMIC signatures in the lymphoma mutational
    catalogues, no cutoff for the LCD (Linear Combination Decomposition)."

## ----exposures_barplot_abs_cutoffZero, fig.width=6, fig.height=5, fig.cap=cap----
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffZero_LCDlist$exposures,
  in_signatures_ind_df = CosmicValid_cutoffZero_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)

## ----caption_barplot_4, echo=FALSE--------------------------------------------
cap="Relative exposures of the COSMIC signatures in the lymphoma mutational
    catalogues, no cutoff for the LCD (Linear Combination Decomposition)."

## ----exposures_barplot_rel_cutoffZero, fig.width=6, fig.height=5, fig.cap=cap----
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffZero_LCDlist$norm_exposures,
  in_signatures_ind_df = CosmicValid_cutoffZero_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)

## ----definition_of_cutoff-----------------------------------------------------
my_cutoff <- 0.06

## ----LCD_complex_cutoff_cutoffGen, warning=FALSE------------------------------
general_cutoff_vector <- rep(my_cutoff,dim(current_sig_df)[2])
CosmicValid_cutoffGen_LCDlist <- LCD_complex_cutoff(
  in_mutation_catalogue_df = lymphomaNature2013_mutCat_df,
  in_signatures_df = current_sig_df,
  in_cutoff_vector = general_cutoff_vector,
  in_sig_ind_df = current_sigInd_df)

## -----------------------------------------------------------------------------
kable(CosmicValid_cutoffGen_LCDlist$out_sig_ind_df, row.names=FALSE,
      caption=paste0("Signatures with cohort-wide exposures > ",my_cutoff))

## ----caption_barplot_5, echo=FALSE--------------------------------------------
cap="Absoute exposures of the COSMIC signatures in the lymphoma mutational
    catalogues, cutoff of 6% for the LCD (Linear Combination Decomposition)."

## ----exposures_barplot_abs_cutoffGen, fig.width=6, fig.height=4, fig.cap=cap----
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffGen_LCDlist$exposures,
  in_signatures_ind_df = CosmicValid_cutoffGen_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)

## ----caption_barplot_6, echo=FALSE--------------------------------------------
cap="Relative exposures of the COSMIC signatures in the lymphoma mutational
    catalogues, cutoff of 6% for the LCD (Linear Combination Decomposition)."

## ----exposures_barplot_rel_cutoffGen, fig.width=6, fig.height=4, fig.cap=cap----
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffGen_LCDlist$norm_exposures,
  in_signatures_ind_df = CosmicValid_cutoffGen_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)

## ----make_specific_cutoff_vector----------------------------------------------
specific_cutoff_vector <- general_cutoff_vector
specific_cutoff_vector[c(1,5)] <- 0
specific_cutoff_vector

## ----LCD_complex_cutoff_cutoffSpec, warning=FALSE-----------------------------
CosmicValid_cutoffSpec_LCDlist <- LCD_complex_cutoff(
  in_mutation_catalogue_df = lymphomaNature2013_mutCat_df,
  in_signatures_df = current_sig_df,
  in_cutoff_vector = specific_cutoff_vector,
  in_sig_ind_df = current_sigInd_df)

## ----caption_barplot_7, echo=FALSE--------------------------------------------
cap="Absoute exposures of the COSMIC signatures in the lymphoma mutational
    catalogues, cutoff of 6% for the LCD (Linear Combination Decomposition)."

## ----exposures_barplot_abs_cutoffSpec, fig.width=6, fig.height=4, fig.cap=cap----
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffSpec_LCDlist$exposures,
  in_signatures_ind_df = CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)

## ----caption_barplot_8, echo=FALSE--------------------------------------------
cap="Relative exposures of the COSMIC signatures in the lymphoma mutational
    catalogues, cutoff of 6% for the LCD (Linear Combination Decomposition)."

## ----exposures_barplot_rel_cutoffSpec, fig.width=6, fig.height=4, fig.cap=cap----
exposures_barplot(
  in_exposures_df = CosmicValid_cutoffSpec_LCDlist$norm_exposures,
  in_signatures_ind_df = CosmicValid_cutoffSpec_LCDlist$out_sig_ind_df,
  in_subgroups_df = COSMIC_subgroups_df)

## ----caption_heatmap, echo=FALSE----------------------------------------------
cap="Clustering of Samples and Signatures based on the relative exposures of
    the COSMIC signatures in the lymphoma mutational catalogues."

## ----apply_heatmap_exposures, fig.width=6, fig.height=4, fig.cap=cap----------
complex_heatmap_exposures(CosmicValid_cutoffGen_LCDlist$norm_exposures,
                          COSMIC_subgroups_df,
                          CosmicValid_cutoffGen_LCDlist$out_sig_ind_df,
                          in_data_type="norm exposures",
                          in_subgroup_colour_column="col",
                          in_method="manhattan",
                          in_subgroup_column="subgroup")

## ----caption_hclust, echo=FALSE-----------------------------------------------
cap="Clustering of the Samples based on the relative exposures of the
    COSMIC signatures in the lymphoma mutational catalogues."

## ----apply_hclust_exposures, fig.width=6, fig.height=3, fig.cap=cap-----------
hclust_list <- 
  hclust_exposures(CosmicValid_cutoffGen_LCDlist$norm_exposures,
                   COSMIC_subgroups_df,
                   in_method="manhattan",
                   in_subgroup_column="subgroup")

## ----caption_heatmap_2, echo=FALSE--------------------------------------------
cap=paste0("PIDs labelled by the clusters extracted from the
            signature analysis.")

## ----cluster_PIDs_sig_info, fig.width=6, fig.height=4, fig.cap=cap------------
cluster_vector <- cutree(hclust_list$hclust,k=4)
COSMIC_subgroups_df$cluster <- cluster_vector
subgroup_colour_vector <- rainbow(length(unique(COSMIC_subgroups_df$cluster)))
COSMIC_subgroups_df$cluster_col <- 
  subgroup_colour_vector[factor(COSMIC_subgroups_df$cluster)]
complex_heatmap_exposures(CosmicValid_cutoffGen_LCDlist$norm_exposures,
                          COSMIC_subgroups_df,
                          CosmicValid_cutoffGen_LCDlist$out_sig_ind_df,
                          in_data_type="norm exposures",
                          in_subgroup_colour_column="cluster_col",
                          in_method="manhattan",
                          in_subgroup_column="cluster")

