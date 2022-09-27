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

kablize <- function(tab, digits = 3) {
  kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}

format_num <- function(number, digits = 0, format = "f") {
  formatC(number, digits = digits, format = format, big.mark = ",")
}

# ######################## Read in Raw IDAT Files #######################
folder_raw_dat <- sample(searchIDATprefixes(here("DataRaw/methylation/RawIdat/")))
#Raw Data
methylation_raw <- lapply(folder_raw_dat,readIDATpair)

#Prep QualityMask -> inferinfiniumIchannel -> POOBAH -> Noob -> DyeBiasNL 
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

mvals <- BetaValueToMValue(betas) %>% 
  rownames_to_column(var = "CpG_Site")

betas <- betas %>% 
  rownames_to_column(var = "CpG_Site")

write_tsv(betas, here("DataProcessed/methylation/methylation_betas_final_2022_09_27.txt"))
write_tsv(mvals, here("DataProcessed/methylation/methylation_mvals_final_2022_09_27.txt"))

