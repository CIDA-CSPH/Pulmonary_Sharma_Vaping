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


kablize <- function(tab, digits = 3) {
  kableExtra::kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}


# Read in files ---------------------------------------------------------------

clin_metadata <- read.csv(here("DataProcessed/methylation/clin_metadata_w_ruv.csv")) %>% 
  dplyr::rename("ruv_k1" = X1,
                "ruv_k2" = X2) %>% 
  drop_na(vape_6mo_lab, methylation_id) %>% 
  mutate(vape_status = factor(vape_6mo_lab, 
                              levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"),
                              labels =  c("Not Vaped", "Vaped")),
         recruitment_center = factor(recruitment_center, 
                                     levels = c("Pueblo", "Aurora", "CommCity/Denver")))

mvals <- read_tsv(here("DataProcessed/methylation/methylation_mvals_final_2022_09_27.txt")) %>% 
  column_to_rownames("CpG_Site") 

mvals <- mvals[,colnames(mvals) %in% clin_metadata$sentrix_name]


# Write a function to fit the model ---------------------------------------
mod_fit <- function(y){
  fit <- lm(y ~ age + recruitment_center + sex_lab + vape_status + ruv_k1 + ruv_k2, data = clin_metadata)
  fit.out <-  summary(fit)
  res <-  fit.out$coefficients["vape_statusVaped",]
  return(res)
}

#start a multisession
plan(multisession)

CpG_Results <- future_apply(mvals, 1, mod_fit)

plan(sequential)

CpG_Results_df <- CpG_Results %>% 
  t() %>% 
  as.data.frame()

CpG_Results_df <- CpG_Results_df %>% 
  mutate(fdr = p.adjust(as.numeric(CpG_Results_df$`Pr(>|t|)`), method = 'fdr')) %>% 
  rename("p.value" = `Pr(>|t|)`) %>% 
  rownames_to_column(var = "CpG_Site")

write_csv(CpG_Results_df, here("DataProcessed/methylation/results/results_ruvk1_ruvk2.csv"))

sig_sites <- CpG_Results_df[CpG_Results_df$fdr < 0.05,]

CpG_Results_df %>% 
  ggplot(aes(x = p.value)) +
  geom_histogram(col = "white", bins = 20)+
  labs(title = "P-value distribution",
       y = "Count",
       x = "p-value")


# Try dropping ruv_k2 -----------------------------------------------------
mod_fit_ruv_k1 <- function(y){
  fit <- lm(y ~ age + recruitment_center + sex_lab + vape_status + ruv_k1, data = clin_metadata)
  fit.out <- summary(fit)
  res <- fit.out$coefficients["vape_statusVaped",]
  return(res)
}

plan(multisession)

ruv_k1_res <- future_apply(mvals, 1, mod_fit_ruv_k1)

ruv_k1_res_df <- ruv_k1_res %>% 
  t() %>% 
  as.data.frame()


ruv_k1_res_df <- ruv_k1_res_df %>% 
  mutate(fdr = p.adjust(as.numeric(ruv_k1_res_df$`Pr(>|t|)`), method = 'fdr')) %>% 
  rename("p.value" = `Pr(>|t|)`)

#write_csv(CpG_Results_df, here("DataProcessed/methylation/results/results_ruvk1.csv"))
