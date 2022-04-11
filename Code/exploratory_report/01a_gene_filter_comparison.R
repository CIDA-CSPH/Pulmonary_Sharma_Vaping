library(tidyverse)
library(here)
library(kableExtra)


#read in raw gene counts
raw_gene_count <- read_tsv(file = here("DataRaw/gene_counts_choo.txt"), col_names = T)

#table to store results
gene_filter_tib <- tibble(Filter = c(1,2,3),
                          Incluison_Criteria = c("At least 25% of the samples have > 0 reads",
                                                   "The range of reads across all samples < 100", 
                                                   ">5 reads in at least 2 samples (Bioconductor)"),
                          Gene_Count_Before = rep(NA, 3),
                          Gene_Count_After = rep(NA, 3),
                          Genes_Removed = rep(NA, 3))



#remove all genes with 0 counts in 75% or more of samples
filter_0_.75_count <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) sum(x == 0) >= 49*.75)

#True if # of cells with 0's >= 75% of samples
#Keep if NOT true
zero_count_.75_filtered <- raw_gene_count[!filter_0_.75_count,]

#remove all genes with range >= 100
filter_range_100 <- zero_count_.75_filtered %>% 
  select(-Feature) %>% 
  apply(.,1, function(x) (max(x) - min(x)) >= 100)
#True if range of read counts across all samples are >= 100
#Keep if TRUE. Want to remove genes with low variation across samples
original_choo_filtered <- zero_count_.75_filtered[filter_range_100,]


#Bioconductor suggested filter (at least two samples with > 5 counts)
biocon_filter <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) length(x[x>5])>=2)
#True if there are at least two samples with >5 reads
#Keep genes if TRUE
bioconductor_filtered <- raw_gene_count[biocon_filter,]

#Assign results to table
gene_filter_tib <- gene_filter_tib %>% 
  dplyr::mutate(Gene_Count_Before = c(nrow(raw_gene_count), nrow(zero_count_.75_filtered), nrow(raw_gene_count)),
                Gene_Count_After = c(nrow(zero_count_.75_filtered), nrow(original_choo_filtered), nrow(bioconductor_filtered)),
                Genes_Removed = c(sum(filter_0_.75_count), sum(!filter_range_100), sum(!biocon_filter)))

filter_compar <- kable(gene_filter_tib, digits = 0, format.args = list(big.mark = ","), 
                              col.names = gsub("_", " ", names(gene_filter_tib)), booktabs = T) %>% 
  kable_styling(latex_options = "striped")

##gene range filter
get_range <- (function(x) (max(x) - min(x)))

#Get the range for each gene
genes_range_after_0_rm <- zero_count_.75_filtered %>% 
  select(-Feature) %>% 
  apply(.,1, get_range)

#Set a range of thresholds
range_to_test <- seq(50,200,10)


#make empty table
range_test_mat <- tibble(range = range_to_test,
                         genes_removed = rep(NA, length(range_to_test)))

#Apply the range of thresholds
gene_gets_removed <- t(sapply(genes_range_after_0_rm, function(x) x < range_to_test))

range_test_mat$genes_removed <- apply(gene_gets_removed, 2, function(x) sum(x))

#Make subset for cutoff used
cutoff_used <- range_test_mat %>% 
  filter(range == 100)

#Show cutoffhow to 
range_cutoff <- range_test_mat %>% 
  ggplot(aes(x = range, y = genes_removed)) +
  geom_point() + 
  geom_point(data = range_test_mat %>% filter(range == 100), 
             aes(color = 'Cutoff Used')) +
  geom_text(data = range_test_mat %>% filter(range == 100), label = "14,645", vjust = -1, color = 'red') +
  labs(x = "Range Threshold", y = "Genes Removed", title = "Genes Removed for Range Filter Cutoff", color = "")

