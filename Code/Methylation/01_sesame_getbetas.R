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
library(BiocParallel)

# ######################## Read in Raw IDAT Files #######################
folder_raw_dat <- searchIDATprefixes(here("DataRaw/methylation/RawIdat/"))

#Raw Data
methylation_raw <- lapply(folder_raw_dat,readIDATpair)

#Prep QualityMask -> inferinfiniumIchannel -> POOBAH -> Noob -> BMIQ -> DyeBiasNL 
#non-experimental masking
methyl_mask <- lapply(methylation_raw, qualityMask)

#infer infinium I channel
methyl_infer <- lapply(methyl_mask, inferInfiniumIChannel)

# Remove low detection p-values and masked probes -------------------------

##get detection pvalues
pvals <- lapply(methyl_infer, pOOBAH, return.pval=T) %>% as.data.frame

##keep pvalues with avg <0.05 across all samples
cpg_keep <- apply(pvals, 1,function (x) {
  sum(x > 0.05) < 0.1*ncol(pvals)
  })

##Keep cpg_keep and unmasked probes
cpg_keep_final <- cpg_keep & !methyl_mask$`205310770093_R03C01`$mask

## Subset the probes for each sample
methyl_db_thin <- lapply(methyl_infer, function (sdf) {
  sdf[cpg_keep_final,]
})


# continue the normalization process --------------------------------------
#noob
methyl_noob <- lapply(methyl_db_thin, noob)

#BMIQ
methyl_bmiq <- lapply(methyl_noob, matchDesign)

#dye bias correction
methyl_norm <- lapply(methyl_bmiq, dyeBiasNL)

#get betas
betas <- do.call(cbind, lapply(methyl_norm, getBetas)) %>% as.data.frame()

mvals_write <- BetaValueToMValue(betas) %>% 
  rownames_to_column(var = "CpG_Site")

betas_write <- betas %>% 
  rownames_to_column(var = "CpG_Site")

#write_tsv(betas_write, here("DataProcessed/methylation/QC/methylation_betas_final_2022_09_27.txt"))
#write_tsv(mvals_write, here("DataProcessed/methylation/QC/methylation_mvals_final_2022_09_27.txt"))


# Autosomal Betas ---------------------------------------------------------
## This function takes your normalized data and gets the betas for only the autosomal probes
get_auto_vals <- function(normal_dat, get_value = "Betas"){
  
  #Get annotation
  autosomes <- sesameData_getAutosomeProbes("EPIC")
  
  cpg_sites <- names(autosomes)
  
  probe_annos <- autosomes %>% as.data.frame(row.names = names(autosomes))
  
  autosomal_probes <- probe_annos[!probe_annos$seqnames %in% c("chrX", "chrY", "chrM"),]
  
  #Filter for only autosomal probes
  methyl_norm_auto <- lapply(normal_dat, function(x) {
    x %>% 
      filter(Probe_ID %in% rownames(autosomal_probes))
  })
  #Noob Only 
  betas_auto <- do.call(cbind, lapply(methyl_norm_auto, getBetas)) %>% as.data.frame()
  
  return(betas_auto)
}

# Get Autosomal betas and write them out to another file ------------------
betas_auto <- get_auto_vals(methyl_norm)

betas_auto_write <- betas_auto %>% 
  rownames_to_column(var = "CpG_Site")

write_tsv(betas_auto_write, here("DataProcessed/methylation/QC/autosomal_betas_final_2022_09_30.txt"))