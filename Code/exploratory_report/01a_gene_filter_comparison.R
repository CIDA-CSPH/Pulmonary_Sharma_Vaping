#load libraries
library(tidyverse)
library(readr)
library(here)
library(kableExtra)

#read in raw gene counts
raw_gene_count <- read_tsv(file = here("DataRaw/gene_counts_choo.txt"), col_names = T)

#table to store results
gene_filter_tab <- matrix(nrow = 4, ncol = 4)

colnames(gene_filter_tab) <- c("Inclusion Criteria", "Read Count Before", "Read Count After", "Reads Removed")

gene_filter_tab[,1] <- c("At least one sample has > 0 reads", "At least 25% of the samples have > 0 reads",
                         "The range of reads across all samples < 100", ">=5 reads in at least 2 samples (Bioconductor)")

gene_filter_tab[,2] <- rep(60651, 4)

#remove all genes with 0 counts
filter_0_count <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) sum(x) == 0)

sum(filter_0_count)

zero_genes_filtered <- raw_gene_count[!filter_0_count,]



#remove all genes with 0 counts in 75% or more of samples
filter_0_.75_count <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) sum(x == 0) >= 49*.75)

sum(filter_0_.75_count)

zero_count_.75_filtered <- raw_gene_count[!filter_0_.75_count,]

#remove all genes with range >= 100
filter_range_100 <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(.,1, function(x) (max(x) - min(x)) >= 100)

sum(filter_range_100)

range_100_filtered <- raw_gene_count[!filter_range_100,]

#choo original code 
#remove counts with range > 100? (Trent)
filter_genes_range <-
  raw_gene_count[ , -1] %>%
  apply(., 1, function(x) { max(x) - min(x) } ) %>%
  `>=`(., 100)
sum(filter_genes_range)

filter_genes_range == filter_range_100

#Bioconductor suggested filter (at least two samples with > 5 counts)
filter <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) length(x[x>5])>=2)

bioconductor_filter <- raw_gene_count[filter,]

gene_count_after <- c(nrow(zero_genes_filtered), nrow(zero_count_.75_filtered), nrow(range_100_filtered), nrow(bioconductor_filter))

gene_filter_tab[,3] <- gene_count_after

gene_filter_tab[,4] <- c(sum(filter_0_count), sum(filter_0_.75_count), sum(filter_range_100), sum(!filter))

kbl(gene_filter_tab, booktabs = T, digits = 3) %>%
  kable_styling(position = "center", latex_options = c("striped", "hold_position"))
