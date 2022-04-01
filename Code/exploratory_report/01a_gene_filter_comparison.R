library(tidyverse)
library(here)


#read in raw gene counts
raw_gene_count <- read_tsv(file = here("DataRaw/gene_counts_choo.txt"), col_names = T)

#table to store results
gene_filter_mat <- matrix(nrow = 4, ncol = 4)

colnames(gene_filter_mat) <- c("Inclusion Criteria", "Read Count Before", "Read Count After", "Reads Removed")

gene_filter_mat[,1] <- c("At least one sample has > 0 reads", "At least 25% of the samples have > 0 reads",
                         "The range of reads across all samples < 100", ">5 reads in at least 2 samples (Bioconductor)")

#remove all genes with 0 counts
filter_0_count <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) sum(x) == 0)


zero_genes_filtered <- raw_gene_count[!filter_0_count,]

#remove all genes with 0 counts in 75% or more of samples
filter_0_.75_count <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) sum(x == 0) >= 49*.75)


zero_count_.75_filtered <- raw_gene_count[!filter_0_.75_count,]

#remove all genes with range >= 100
filter_range_100 <- zero_count_.75_filtered %>% 
  select(-Feature) %>% 
  apply(.,1, function(x) (max(x) - min(x)) >= 100)

range_100_filtered <- zero_count_.75_filtered[!filter_range_100,]


#Bioconductor suggested filter (at least two samples with > 5 counts)
filter <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) length(x[x>5])>=2)

bioconductor_filter <- raw_gene_count[filter,]

gene_count_after <- c(nrow(zero_genes_filtered), nrow(zero_count_.75_filtered), nrow(range_100_filtered), nrow(bioconductor_filter))

gene_filter_mat[,2] <- c(nrow(raw_gene_count), nrow(raw_gene_count), nrow(zero_count_.75_filtered), nrow(raw_gene_count))

gene_filter_mat[,3] <- gene_count_after

gene_filter_mat[,4] <- c(sum(filter_0_count), sum(filter_0_.75_count), sum(filter_range_100), sum(!filter))

gene_filter_tib <- as_tibble(gene_filter_mat) %>% 
  filter(gene_filter_mat[,3] != 46356) %>% 
  mutate(Filter = c(1,2,3), "Analysis Used" = c("Previous", "Previous", "Current")) %>% 
  select(Filter, "Analysis Used", everything())

filter_compar <- kbl(gene_filter_tib, booktabs = T, digits = 3) %>%
  kable_styling(position = "center", latex_options = c("striped", "hold_position"))

##gene range filter
range_filter <- (function(x) (max(x) - min(x)))

#Get the range for each gene
filter_range_test <- zero_count_.75_filtered %>% 
  select(-Feature) %>% 
  apply(.,1, function(x) (max(x) - min(x)))

#Set a range of thresholds
range_to_test <- seq(50,200,10)

#Apply the range of thresholds
genes_removed <- t(sapply(filter_range_test, function(x) x <= range_to_test)) %>% 
  apply(., 2, function(x) sum(x))


#make empty table
range_test_mat <- tibble(range = range_to_test,
                         genes_removed = genes_removed)

#Make subset for cutoff used
cutoff_used <- range_test_mat %>% 
  filter(range == 100)

#Show cutoffhow to 
range_cutoff <- range_test_mat %>% 
  ggplot(aes(x = range, y = genes_removed)) +
  geom_point() + 
  geom_point(data = range_test_mat %>% filter(range == 100), 
             aes(color = 'Cutoff Used')) +
  geom_text(data = range_test_mat %>% filter(range == 100), label = "14678", vjust = -1, color = 'red') +
  labs(x = "Range Threshold", y = "Genes Removed", title = "Genes Removed for Range Filter Cutoff", color = "Legend")

