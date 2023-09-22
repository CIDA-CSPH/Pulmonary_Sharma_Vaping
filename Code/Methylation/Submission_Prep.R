## ---------------------------
## Script name: Formatting data for paper submission
##
## Purpose of script:
##
## Author: Trent Hawkins
##
## Date Created: 2023-05-03
## ---------------------------

# Methylation -------------------------------------------------------------

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

metadata <- read_csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv"))

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

pvals_join <- pvals %>% rownames_to_column(var = "CpG_Site")

colnames(pvals_join) <- gsub("X", "", colnames(pvals_join))

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

betas_write <- betas %>% 
  rownames_to_column(var = "CpG_Site")


# Map sentrix name to Sample ID -------------------------------------------
all(colnames(betas_write) == colnames(pvals_join))

sentrix <- colnames(betas) %>% as.data.frame
colnames(sentrix) <- "sentrix_name"

sentrix_map <- left_join(sentrix, metadata %>% select(rna_id, sentrix_name), by = "sentrix_name")

## Change names in betas and pvals
colnames(betas_write) <- c("CpG_Site", sentrix_map$rna_id)
colnames(pvals_join) <- c("CpG_Site", sentrix_map$rna_id)
# Join betas and pvals ----------------------------------------------------
geo_sub <- left_join(pvals_join, betas_write, by = "CpG_Site", suffix = c(".beta", ".pval"))

colnames(geo_sub) <- c("CpG_Site", sort(colnames(geo_sub)[-1]))

samples_map <- metadata %>% 
  select(rna_id, sentrix_name) %>% 
  arrange(sentrix_name) %>% 
  drop_na()

new_colnames <- c(rbind(samples_map$rna_id, "Detection Pval"))

colnames(geo_sub) <- c("CpG_Site", new_colnames)

head(geo_sub)

write_tsv(geo_sub, here("DataProcessed/methylation/Submission/betas_detPvals.txt"))


