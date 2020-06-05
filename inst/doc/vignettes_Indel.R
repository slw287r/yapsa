## ----load_style, warning=FALSE, echo=FALSE, message=FALSE, results="hide"-----
library(BiocStyle)

## ----packages, include=FALSE--------------------------------------------------
library(YAPSA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(knitr)
opts_chunk$set(echo=TRUE)
opts_chunk$set(fig.show='asis')

## ----load_PCAWG_sigs----------------------------------------------------------
data(sigs_pcawg)

## ----caption_spectra, echo=FALSE----------------------------------------------
cap <-": Nucleotide exchange spectra of the Indel signatures ID3, associated 
        with tobacco smoking, and ID6, related to deficiencies in homologous 
        recombination repair."

## ----INDEL_sig_example, include=TRUE, fig.width=15, fig.height=6, fig.cap= cap----
plotExchangeSpectra_indel(PCAWG_SP_ID_sigs_df[,c(3,6)])

## ----INDEL_sig_info-----------------------------------------------------------
current_caption <- paste0("Information on Indel mutational signatures.")
if(!exists("repress_tables"))
  kable(PCAWG_SP_ID_sigInd_df, row.names=FALSE, caption=current_caption)

## ----load_GoNL----------------------------------------------------------------
data(GenomeOfNl_raw)
GenomeOfNl_raw <- GenomeOfNl_raw[, c(1,2,4,5)]

## ----load_GoNL_raw------------------------------------------------------------
load_data_new <- FALSE
if(load_data_new){
  data <- data.frame(matrix(ncol = 8, nrow = 0))
  for(index in seq_along(1:22)){
              print(index)
              temp <- tempfile()
              file_path <- paste0("https://molgenis26.target.rug.nl/
                                  downloads/gonl_public/variants/release5/
                                  gonl.chr",
                                  index, ".snps_indels.r5.vcf.gz")
              download.file(file_path, temp)
       
              data <- rbind(data, read.table(gzfile(temp, paste0("gonl.chr",
                                                    index,
                                                    ".snps_indels.r5.vcf")),
                                             header=FALSE, sep="\t", 
                                             stringsAsFactors = FALSE))
              data <- data[grep("INDEL", data$V8),]
              unlink(temp)
            
  }
  colnames(data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                      "FILTER","INFO")

  GenomeOfNl_raw <- data[, c(1,2,4,5)]
}

## ----shows_top_of_df, echo=FALSE----------------------------------------------
kable(head(GenomeOfNl_raw), caption="Head of VCF file 
      containing the GoNL INDEL data")

## ----randomize_data_set, warning= FALSE---------------------------------------
seed = 2
set.seed(seed)
number_of_indels <- sample(c(30:70), 15,  replace = TRUE)

index=0
seed=3
set.seed(seed)
vcf_like_indel_lists <- lapply(number_of_indels, function(size){
  df_per_PID <- GenomeOfNl_raw[sample(nrow(GenomeOfNl_raw), size, 
                                      replace = FALSE), ]
  index <<- index+1
  df_per_PID$PID <- rep(paste0("PID_", index), length(size))
  df_PIDs <- df_per_PID[order(df_per_PID$CHROM),]
  return(df_PIDs)
  })

vcf_like_indel_df <- do.call(rbind.data.frame, vcf_like_indel_lists)
kable(head(vcf_like_indel_df), caption="Head of the vcf_like_df
      containing the subsampled GoNL Indel data")

## ----create_mutational_catalog, warning=FALSE---------------------------------
vcf_like_indel_trans_df <- translate_to_hg19(vcf_like_indel_df,"CHROM")
mutational_cataloge_indel_df <- create_indel_mutation_catalogue_from_df(
  in_dat = vcf_like_indel_trans_df,
  in_signature_df = PCAWG_SP_ID_sigs_df)

kable(head(mutational_cataloge_indel_df[,1:5]))

## ----load_cutoffs, warning=FALSE----------------------------------------------
data(cutoffs_pcawg)

## ----LCD_decompostion, warning=FALSE------------------------------------------
current_catalogue_df <- mutational_cataloge_indel_df 
current_sig_df <- PCAWG_SP_ID_sigs_df
current_cutoff_pid_vector <- cutoffPCAWG_ID_WGS_Pid_df[3,]
current_sigInd_df <- PCAWG_SP_ID_sigInd_df

current_LCDlistsList <- LCD_complex_cutoff_combined(
  current_catalogue_df,
  current_sig_df,
  in_cutoff_vector = current_cutoff_pid_vector,
  in_filename = NULL,
  in_method = "abs",
  in_sig_ind_df = current_sigInd_df)

current_consensus_LCDlist <- current_LCDlistsList$consensus
if(!exists("repress_tables")) 
  as.character(current_consensus_LCDlist$out_sig_ind_df$sig)

## ----caption_exposure, echo=FALSE---------------------------------------------
cap <- ":Exposures to Indel mutational signatures in the artificial data created 
        by sampling GoNL variants. Exposures were obtained from
        a decomposition with PCAWG Indel signatures as well as their signature
        specific-cutoffs (cutoffPCAWG_ID_WGS_Pid_df)."

## ----plot_exposure, echo=TRUE, warning=FALSE, fig.width=15, fig.height=6, fig.cap= cap----
exposures_barplot(current_LCDlistsList$perPID$exposures,
                  current_LCDlistsList$perPID$out_sig_ind_df)

## ----caption_CI, echo=FALSE---------------------------------------------------
cap <- "Confidence interval calculation for exposures to Indel mutational 
        signatures"

## ----CI, echo=TRUE, warning=FALSE, fig.width=17, fig.height=15, fig.cap=cap----
confidence_intervals_ID <- confidence_indel_only_calulation(
  in_current_indel_df = current_catalogue_df)
plot(confidence_intervals_ID$p_complete_PCAWG_ID)

