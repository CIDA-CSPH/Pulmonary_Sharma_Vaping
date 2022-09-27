####################Load Libraries########################
library(tidyverse)
library(here)
library(kableExtra)
library(sesame)
library(readxl)
library(parallel)
library(randomForest)
library(janitor)
library(minfi)
library(data.table)
library(ggpubr)

kablize <- function(tab, digits = 3) {
  kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}

format_num <- function(number, digits = 0, format = "f") {
  formatC(number, digits = digits, format = format, big.mark = ",")
}

#Read in Methylation QC File
methylation_qc <- read_csv(here("DataProcessed/methylation/methylation_qc_metrics.csv"))

methylation_qc <- left_join(methylation_qc, metadata_sex %>% select(rna_id, methylation_id, sid, sentrix_name, vape_6mo_lab), by = "sentrix_name")

# ######################## Read in Raw IDAT Files #######################
set.seed(404)
folder_raw_dat <- sample(searchIDATprefixes(here("DataRaw/methylation/RawIdat/")),3)
#Raw Data
methylation_raw <- lapply(folder_raw_dat,readIDATpair)

# #Prep QualityMask -> inferinfiniumIchannel -> POOBAH -> Noob -> DyeBiasNL 
# methyl_mask <- lapply(methylation_raw, qualityMask)
# 
# methyl_infer <- lapply(methyl_mask, inferInfiniumIChannel)
# 
# methyl_pdet <- lapply(methyl_infer, pOOBAH)
# 
# methyl_noob <- lapply(methyl_pdet, noob)
# 
# methyl_bmiq <- lapply(methyl_noob, matchDesign)
# 
# methyl_noob_dbcorr <- lapply(methyl_noob, dyeBiasNL)
# 
# methyl_bmiq_dbcorr <- lapply(methyl_bmiq, dyeBiasNL)

#Prep QualityMask -> inferinfiniumIchannel -> POOBAH -> Noob -> DyeBiasNL 
methyl_mask <- lapply(methylation_raw, qualityMask)

methyl_infer <- lapply(methyl_mask, inferInfiniumIChannel)

methyl_db <- lapply(methyl_infer, dyeBiasNL)

methyl_pdet <- lapply(methyl_db, pOOBAH, return.pval=T)

methyl_noob <- lapply(methyl_pdet, noob)

methyl_bmiq <- lapply(methyl_noob, matchDesign)

methyl_noob_dbcorr <- lapply(methyl_noob, dyeBiasNL)

methyl_bmiq_dbcorr <- lapply(methyl_bmiq, dyeBiasNL)


# Plot RedGRNQQ -----------------------------------------------------------

# Make redgrnqq plot matrix 
plot_redgrnqq <- function(raw_dat, noob_dat, bmiq_dat){
  
  par(mfrow=c(3,3), mar=c(3,3,2,1))
  
  #Make Raw Plots
  for (h in rand_sdf){
    sesameQC_plotRedGrnQQ(raw_dat[[h]], main= paste0(h, " (Raw)"))
    
    
  }
  #Make noob Plots
  for (i in rand_sdf){
    sesameQC_plotRedGrnQQ(noob_dat[[i]],
                          main= paste0(i, " (noob)"))
  }
  
  #Make noob + BMIQ Plots
  for (j in rand_sdf){
    sesameQC_plotRedGrnQQ(bmiq_dat[[j]],
                          main= paste0(j, " (noob + BMIQ)"))
  }
}

plot_redgrnqq(raw_dat = methylation_raw,
              noob_dat = methyl_noob_dbcorr,
              bmiq_dat = methyl_bmiq_dbcorr)


# Look at betas dist'n ----------------------------------------------------

# Make redgrnqq plot matrix 
plot_betas_compare <- function(noob_dat, bmiq_dat){
  
  par(mfrow=c(2,3))
  
  #Make noob Plots
  for (i in rand_sdf){
    sesameQC_plotBetaByDesign(noob_dat[[i]],
                              main= paste0(i, " (noob)"))
  }
  
  #Make noob + BMIQ Plots
  for (j in rand_sdf){
    sesameQC_plotBetaByDesign(bmiq_dat[[j]],
                              main= paste0(j, " (noob + BMIQ)"))
  }
}

plot_betas_compare(noob_dat = methyl_noob_dbcorr,
              bmiq_dat = methyl_bmiq_dbcorr)
