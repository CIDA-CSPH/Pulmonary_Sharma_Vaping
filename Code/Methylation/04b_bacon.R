## ---------------------------
## Script name: Bacon 
##
## Purpose of script: Running Bacon to correct deflated p-values
##
## Author: Trent Hawkins
##
## Date Created: 2022-11-29
## ---------------------------
## set working directory for Mac and PC

library(tidyverse)
library(here)
library(future.apply)
library(bacon)
library(limma)


kablize <- function(tab, digits = 3) {
  kableExtra::kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}


# Read in the results -----------------------------------------------------


full_res <- read_csv(here("DataProcessed/methylation/results/results_ruvk1_ruvk2.csv"))


# run bacon ---------------------------------------------------------------

bc <- bacon(full_res$`t value`)

esetimates.corr <- estimates(bc)
bc_inflation <- inflation(bc)
bc_bias <- bias(bc)
pvals.corr <- pval(bc)

full_res$pval.bacon <- as.vector(pvals.corr)
full_res$fdr.bacon <- p.adjust(full_res$pval.bacon, method = "fdr")

sum(full_res$fdr.bacon < 0.05)

sig_results <- full_res[full_res$fdr.bacon < 0.05, ]

hist(full_res$pval.bacon,
     main = "P-Value Distribution Using Bacon (Full Model)",
     xlab = "p-value")
