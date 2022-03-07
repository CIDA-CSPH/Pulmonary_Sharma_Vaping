# run script 10 and 11 then 

# 60651 x 49
# convert the IDs given to the core to PIDs
metadata_final <-
  metadata_naive %>%
  filter(SID != "102") %>%
  mutate(AgeScaled = scale(Age))

final_counts <-
  counts_mat_naive %>%
  .[ , colnames(.) %in% metadata_final$NewID]
library(tidyverse)
library(edgeR)

d <- DGEList(final_counts)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
disp <- d$tagwise.dispersion


mu <- apply(final_counts, 1, mean)

library(ssizeRNA)
powercalc <-
  ssizeRNA_vary(
    nGenes = nrow(final_counts), pi0 = 0.8, m = 200, mu = mu,
    disp = disp, fc = c(1.25),
    fdr = 0.05,
    power = 0.8, maxN = 200)

50*0.9
powercalc$power

