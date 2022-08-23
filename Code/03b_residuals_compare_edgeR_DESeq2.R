#load libraries
library(tidyverse)
library(readr)
library(janitor)
library(RUVSeq)
library(here)
library(DESeq2)
library(Hmisc)
library(reshape2)

#get residuals from edgeR Method

edgeR_Resid <- read_csv(here("DataProcessed/rna_seq/ruv/first_pass_residuals_edgeR.csv"))

DESeq_Resid <- read_csv(here("DataProcessed/rna_seq/ruv/first_pass_residuals_DESeq2.csv"))

colnames(edgeR_Resid) <- colnames(DESeq_Resid)

#subtract the residuals between the two
resid_diff <- edgeR_Resid - DESeq_Resid

#Plot them
set.seed(404)
resid_diff_randplot <- resid_diff[,sample(ncol(resid_diff), size = 4), drop = T]

resid_diff_melt <- melt(resid_diff_randplot)

resid_diff_hist <- resid_diff_melt %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~variable) + 
  labs(title = 'Difference in residuals (edgeR - DESeq2)')

