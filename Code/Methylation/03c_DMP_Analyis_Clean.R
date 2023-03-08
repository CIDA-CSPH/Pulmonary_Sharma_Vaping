## ---------------------------
## Script name: Clean DMP Analysis
##
## Purpose of script: To provide clean, understandable code for the results reported in "DMP_Results".
##
## Author: Trent Hawkins
##
## Date Created: 2022-12-01
## ---------------------------


## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(here)
require(bacon)
require(future.apply)

# Read in files ---------------------------------------------------------------

# Read in the  clinical metadata and clean up factors for analysis
clin_metadata <- read.csv(here("DataProcessed/methylation/QC/clin_metadata_w_ruv.csv")) %>% 
  dplyr::rename("ruv_k1" = X1,
                "ruv_k2" = X2) %>% 
  drop_na(vape_6mo_lab, methylation_id) %>% 
  mutate(vape_status = factor(vape_6mo_lab, 
                              levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"),
                              labels =  c("Not Vaped", "Vaped")),
         recruitment_center = factor(recruitment_center, 
                                     levels = c("Pueblo", "Aurora", "CommCity/Denver")),
         sex_lab = factor(sex_lab, levels = c('Female', 'Male')))

#Read in the m-values
mvals <- read_tsv(here("DataProcessed/methylation/methylation_mvals_final_2022_09_27.txt")) %>% 
  column_to_rownames("CpG_Site") 

#Subset the values to include only the samples with complete methylation data and vape status
mvals <- mvals[,colnames(mvals) %in% clin_metadata$sentrix_name]

# Write a function to fit the model ---------------------------------------
full_modfit <- function(y){
  fit <- lm(y ~ age + sex_lab + recruitment_center + vape_status + ruv_k1 + ruv_k2, data = clin_metadata)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

no_center_modfit <- function(y){
  fit <- lm(y ~ age + sex_lab + vape_status + ruv_k1 + ruv_k2, data = clin_metadata)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

plan(multisession)

full_results <- future_apply(mvals, 1, full_modfit)

noCenter_results <- future_apply(mvals, 1, no_center_modfit)


plan(sequential)

# Format the results ---------------------------------------------------
format_res <- function(res) {
  res_df <- res %>% 
    t() %>% 
    as.data.frame()
  
  res_df <- res_df %>% 
    mutate(fdr = p.adjust(as.numeric(res_df$`Pr(>|t|)`), method = 'fdr')) %>% 
    rename("p.value" = `Pr(>|t|)`) %>% 
    rownames_to_column(var = "CpG_Site")
}

full_results <- format_res(full_results)

noCenter_results <- format_res(noCenter_results)

# Bacon for full results --------------------------------------------------

# Running bacon
bc <- bacon(full_results$`t value`)

## Extract p-values from bacon
pvals.corr <- pval(bc)

## Add back to full results 
full_results$pval.bacon <- as.vector(pvals.corr)
full_results$fdr.bacon <- p.adjust(full_results$pval.bacon, method = "fdr")

estimates(bc)
# Write out the results ---------------------------------------------------

#write_csv(full_results, here("DataProcessed/methylation/EWAS/results/full_res_bacon_2022_12_01.csv"))

#write_csv(noCenter_results, here("DDataProcessed/methylation/EWAS/results/noCenter_res_2022_12_01.csv"))
