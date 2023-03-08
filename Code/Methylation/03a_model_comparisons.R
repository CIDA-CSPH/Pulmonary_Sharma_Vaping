## ---------------------------
## Script name: Methylation Results Modeling
##
## Purpose of script:
##
## Author: Trent Hawkins
##
## Date Created: 2022-10-31
## ---------------------------


## load up the packages we will need:  (uncomment as required)

library(tidyverse)
library(here)
library(future.apply)
library(bacon)
library(limma)


kablize <- function(tab, digits = 3) {
  kableExtra::kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}


# Read in files ---------------------------------------------------------------

clin_metadata <- read.csv(here("DataProcessed/methylation/QC/clin_metadata_w_ruv.csv")) %>% 
  dplyr::rename("ruv_k1" = X1,
                "ruv_k2" = X2) %>% 
  drop_na(vape_6mo_lab, methylation_id) %>% 
  mutate(vape_status = factor(vape_6mo_lab, 
                              levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"),
                              labels =  c("Not Vaped", "Vaped")),
         recruitment_center = factor(recruitment_center, 
                                     levels = c("Pueblo", "Aurora", "CommCity/Denver")))

clin_metadata_ruvr <- read.csv(here("DataProcessed/methylation/clin_metadata_w_RUVr.csv")) %>% 
  mutate(vape_status = factor(vape_6mo_lab, 
                              levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"),
                              labels =  c("Not Vaped", "Vaped")),
         recruitment_center = factor(recruitment_center, 
                                     levels = c("Pueblo", "Aurora", "CommCity/Denver")))
  

# Read in Mvalues ---------------------------------------------------------
mvals <- read_tsv(here("DataProcessed/methylation/methylation_mvals_final_2022_09_27.txt")) %>% 
  column_to_rownames("CpG_Site") 

mvals <- mvals[,colnames(mvals) %in% clin_metadata$sentrix_name]

# Write a function to fit the model ---------------------------------------
## Full Model
full_mod_fit <- function(y, dat){
  fit <- lm(y ~ vape_status + age + sex_lab + recruitment_center + ruv_k1 + ruv_k2, data = dat)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

## Reduced Model
reduced_mod_fit <- function(y){
  fit <- lm(y ~ vape_status + age + sex_lab, data = clin_metadata)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

## Reduced + Center
center_mod_fit <- function(y){
  fit <- lm(y ~ vape_status + age + sex_lab + recruitment_center, data = clin_metadata)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

## Reduced + Center + K1
center_ruvk1_mod_fit <- function(y){
  fit <- lm(y ~ vape_status + age + sex_lab + recruitment_center + ruv_k1, data = clin_metadata)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

## Reduced + Center + K2
center_ruvk2_mod_fit <- function(y){
  fit <- lm(y ~ vape_status + age + sex_lab + recruitment_center + ruv_k2, data = clin_metadata)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

## RUV only 
ruv_mod_fit <- function(y){
  fit <- lm(y ~ vape_status + age + sex_lab + ruv_k1 + ruv_k2, data = clin_metadata)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

## Formatting results for plotting
format_res <- function(res) {
  res_df <- res %>% 
    t() %>% 
    as.data.frame()
  
  res_df <- res_df %>% 
    mutate(fdr = p.adjust(as.numeric(res_df$`Pr(>|t|)`), method = 'fdr')) %>% 
    rename("p.value" = `Pr(>|t|)`) %>% 
    rownames_to_column(var = "CpG_Site")
}

## p-value histograms
pvalHist <- function(res, name){
 hist(res$p.value, breaks = seq(0, 1, 0.05),
      main = name,
      xlab = "p-value")
}

## pval histogram full function
getPvalHist <- function(mod, name) {
  pvalHist(format_res(future_apply(mvals, 1, mod)), name)
}

# Run models with RUV factors from MissMethyl and get p-value hist --------
par(mfrow = c(2, 3))
plan(multisession)
getPvalHist(reduced_mod_fit, name = "No extra Covariates")
getPvalHist(ruv_mod_fit, name = "RUV_K1 + RUV_K2")
getPvalHist(center_mod_fit, name = "Recruitment Center")
getPvalHist(center_ruvk1_mod_fit, name = "Recruitment Center + RUV_K1")
getPvalHist(center_ruvk2_mod_fit, name = "Recruitment Center + RUV_K2")
getPvalHist(full_mod_fit, name = "Full Model")



# Run model with RUVr  ----------------------------------------------------
plan(multisession)
RUVr_Results <- future_apply(mvals, 1, function(y) {full_mod_fit(y = y, dat = clin_metadata)})

## Format results
RUVr_Results <- format_res(RUVr_Results)

## P-value histogram
hist(RUVr_Results$p.value)

## Write out the results for report
#write_csv(RUVr_Results, here("DataProcessed/methylation/results/results_RUVr_k2.csv"))
