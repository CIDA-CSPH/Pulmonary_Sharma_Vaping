## ---------------------------
## Script name: DESEQ Results Formatting
##
## Purpose of script: Fomrat results for Sarah 
##
## Author: Trent Hawkins
##
## Date Created: 2022-12-01
## ---------------------------
## set working directory for Mac and PC

require(tidyverse)
require(here)

results <- read_csv(here("DataProcessed/rna_seq/differential_expression/full_analysis/full_vape_res_2022_10_13.csv"))

results <- results %>%
  mutate(Signif = if_else(padj < 0.05, "*", "")) %>% 
  rename("ENSG" = ensg,
         "Symbol" = gene_name,
         "Log2FoldChange" = log2FoldChange,
         "p-value" = pvalue,
         "FDR" = padj) %>% 
  select(-c(gene_type, baseMean, lfcSE, stat))

write_csv(results, here("DataProcessed/rna_seq/differential_expression/full_analysis/vape_results_reporting_2022_12_01.csv"))
