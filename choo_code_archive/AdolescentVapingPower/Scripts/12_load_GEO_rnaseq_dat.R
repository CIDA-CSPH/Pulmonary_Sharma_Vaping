library(tidyverse)
library(GEOquery)

library(data.table)


geo_obj <- GEOquery::getGEO('GSE152004')
geo_metadata <-
  geo_obj$GSE152004_series_matrix.txt.gz %>% pData

filter_subj_HC <-
  geo_metadata$`asthma status:ch1` == "healthy control"

sum(filter_subj_HC)
rawcounts <-
  fread("Data/GSE152004_695_raw_counts.txt")
totcounts <-
  rawcounts[ , -1] %>% apply(., 2, sum) %>% `/`(., 1e6)
summary(totcounts)


counts_HC <-
  rawcounts %>%
  select(., which(filter_subj_HC) + 1) # plus one because first row = gene names

dim(counts_HC)
totcounts_HC <-
  totcounts[filter_subj_HC]

mean_gene_cpm_HC <-
  sapply(1:nrow(counts_HC),
         function(i) {
           cpm <- as.numeric(counts_HC[i, -1]) / totcounts_HC
           mean(cpm)
            })
filter_gene_cpm <-
  mean_gene_cpm > 0.5
