#load libraries
library(tidyverse)
library(readr)
library(RUVSeq)
library(here)
library(edgeR)

#read in raw gene counts
raw_gene_count <- read_tsv(file = here("DataRaw/gene_counts_choo.txt"), col_names = T)

## Filtering
# #remove 0 counts (sum over rows)
# #find the genes with 0 count
# raw_gene_count <- raw_gene_count %>% 
#   rowwise() %>% 
#   mutate(sum_to_0 = sum(c_across(-Feature)) == 0)
# 
# #check how many genes that removes
# raw_gene_count %>% 
#   group_by(sum_to_0) %>% 
#   summarise(N = n())

#check Choo's original filter (0 count for more than 75% of samples)
raw_gene_count <- raw_gene_count %>% 
  rowwise() %>% 
  mutate(sum_to_0_75 = sum(c_across(-Feature) == 0) >= 49*0.75)

#check how many genes that removes
raw_gene_count %>% 
  group_by(sum_to_0_75) %>% 
  summarise(N = n())

#filter them out
raw_gene_count <- raw_gene_count %>% 
  filter(sum_to_0_75 == F) %>% 
  select(-sum_to_0_75)

#Check Choo's range filter
raw_gene_count <- raw_gene_count %>% 
  rowwise() %>% 
  mutate(range_above_100 = (max(c_across(-Feature) - min(c_across(-Feature)))) >= 100)

#Check how many genes that removes
raw_gene_count %>% 
  group_by(range_above_100) %>% 
  summarise(N = n())

#Remove the genes
raw_gene_count <- raw_gene_count %>% 
  filter(range_above_100 == T) %>% 
  select(-range_above_100)

##Prepare for RUV
#get residuals matrix



ruv_results_k3 <- RUVr(x = raw_gene_count,
                    k = 3,)