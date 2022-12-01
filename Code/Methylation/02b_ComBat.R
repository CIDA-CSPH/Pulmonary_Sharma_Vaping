## ---------------------------
## Script name: ComBat
##
## Purpose of script: Remove "site" effect from data
##
## Author: Trent Hawkins
##
## Date Created: 2022-11-18
## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(here)
require(sva)
require(BCconf)
require(future.apply)


# Read in Clinical Data ---------------------------------------------------
clin_metadata <- read_csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv")) %>% 
  drop_na(methylation_id, vape_6mo_lab)

# Read in Normalized M-values ---------------------------------------------
mvals <- read_tsv(here("DataProcessed/methylation/methylation_mvals_final_2022_09_27.txt")) %>% 
  column_to_rownames("CpG_Site") 

## Drop Individual missing Vape status
mvals <- mvals[,colnames(mvals) %in% clin_metadata$sentrix_name]


# Order the clin_metadata to match the mvals ------------------------------
clin_metadata <- clin_metadata[match(names(mvals), clin_metadata$sentrix_name),]

## Check that samples are in order
all(clin_metadata$sentrix_name == names(mvals))

# Specify The Model Matrix ------------------------------------------------
modcombat <- model.matrix(~ vape_6mo_lab + age + sex_lab, data = clin_metadata) #Treating recruitment center as batch effect here, so left out of model

center <- clin_metadata$recruitment_center
# Run ComBat --------------------------------------------------------------
combat_mvals <- ComBat(dat = mvals, batch = center, prior.plots = T)

combat_mvals <- combat_mvals %>% as.data.frame()

#Write out the results
#write_csv(combat_mvals, here("DataProcessed/methylation/combat_mvals.csv"))


# Latent Factor Analaysis with BCConf -------------------------------------

bcconfFactors <- Correction(Y = combat_mvals %>% as.matrix(), X = modcombat, ind.cov = c(2), r.confound = 2, method = "huber")

latentFactors <- bcconfFactors$C %>% as.data.frame()
names(latentFactors) <- c("BC_1", "BC_2")

clin_metadata_BCfactors <- cbind(clin_metadata, latentFactors)

#write_csv(clin_metadata_BCfactors, here("DataProcessed/clinical_metadata/clin_metadata_combat_bcfactors.csv"))


# Get model to fit --------------------------------------------------------
clin_metadata_BCfactors <- clin_metadata_BCfactors %>% 
  mutate(vape_status = if_else(vape_6mo_lab == "Did Not Vape in Last 6 Months", "Not Vaped", "Vaped"))

mod_fit <- function(y){
  fit <- lm(y ~ age + sex_lab + vape_status + BC_1 + BC_2, data = clin_metadata_BCfactors)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

plan(multisession)
full_res <- future_apply(combat_mvals, 1, mod_fit)

format_res <- function(res) {
  res_df <- res %>% 
    t() %>% 
    as.data.frame()
  
  res_df <- res_df %>% 
    mutate(fdr = p.adjust(as.numeric(res_df$`Pr(>|t|)`), method = 'fdr')) %>% 
    rename("p.value" = `Pr(>|t|)`) %>% 
    rownames_to_column(var = "CpG_Site")
}

full_res <- format_res(full_res)

hist(full_res$p.value, 
     main = "ComBat -> BCConf (k = 2 factors)",
     xlab = "p-value")

ruvr_res <- read_csv(here("DataProcessed/methylation/results/results_RUVr_k2.csv"))

hist(ruvr_res$p.value,
     main = "Full Model With k = 2 RUVr Factors",
     xlab = "p-value")
