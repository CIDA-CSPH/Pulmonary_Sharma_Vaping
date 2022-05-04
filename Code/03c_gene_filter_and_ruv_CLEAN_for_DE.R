#load libraries
library(tidyverse)
library(readr)
library(RUVSeq)
library(here)
library(edgeR)
library(janitor)


#read in raw gene counts
raw_gene_count <- read_tsv(file = here("DataRaw/gene_counts_choo.txt"), col_names = T)


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
filtered_gene_count <- zero_count_.75_filtered[filter_range_100,]

#convert gene count matrix
filtered_gene_count <- as.data.frame(filtered_gene_count)
rownames(filtered_gene_count) <- filtered_gene_count[,1]
filtered_gene_count <- filtered_gene_count[,-1]

#create a list of genes
genes <- rownames(filtered_gene_count)

#Load metadata
id_relate <- read_tsv(file = here("DataRaw/20201216_coreID_to_PID.txt"), col_names = T) %>% clean_names()
metadata_unjoined <- read_csv(file = here("DataProcessed/metadata_cleaning/table1_clean_data_2022_03_30.csv"))


#Join metadata
metadata_joined <- id_relate %>% 
  mutate(new_id = str_pad(new_id, 2, "left", "0") %>% paste0("Sample", .)) %>% 
  left_join(metadata_unjoined, by = "sid") %>% 
  filter(new_id %in% names(filtered_gene_count)) %>% 
  as.data.frame()

#Filter out subjects missing vape status
metadata_joined <- metadata_joined %>% 
  filter(!is.na(vape_6mo_lab))

#Make sure the count data matches
filtered_gene_count <- filtered_gene_count[,metadata_joined$new_id]


##Prepare for DESeq2
#Set up factors properly
metadata_joined$sex_lab <- factor(metadata_joined$sex_lab, 
                                  levels = c("Female", "Male"))
metadata_joined$vape_6mo_lab <- factor(metadata_joined$vape_6mo_lab, 
                                       levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"),
                                       labels = c("Did_Not_Vape_in_Last_6_Months", "Vaped_in_Last_6_Months"))


#check that row and column names are equal
#metadata_joined$new_id == colnames(filtered_gene_count)

##Prepare for RUV
#make easier to reference
vape_status <- metadata_joined$vape_6mo_lab

male <- metadata_joined$sex_lab

age <- metadata_joined$age

#convert to SEqExpressionSet object 
ruv_prep <- newSeqExpressionSet(as.matrix(filtered_gene_count), 
                                   phenoData = data.frame(vape_status, 
                                                          male, 
                                                          age, 
                                                          row.names = metadata_joined$new_id))


#Read in residuals matrix (See "03_gene_filter_and_ruv_edgeR.R" for more information on calculation)
first_pass_residuls <- as.matrix(read_csv(here("DataProcessed/first_pass_residuals_edgeR.csv")))
rownames(first_pass_residuls) <- genes
colnames(first_pass_residuls) <- metadata_joined$new_id

## Run RUVr with k = 2
ruv_k2 <- RUVr(ruv_prep, cIdx = genes, k = 2, residuals = first_pass_residuls)

#write out needed results
# ruv_k2_counts <- t((counts(ruv_k2[goi, ])+.5)) %>%
#   merge(metadata_joined, ., by="row.names") %>%
#   gather(gene, expression, (ncol(.)-length(goi) + 1):ncol(.))

#write_csv(ruv_k2_counts, file = here("DataProcessed/RUV_k2_Counts_2022_05_01.csv"))

ruv_k2_write <- pData(ruv_k2) %>% 
  mutate(new_id = rownames(.)) %>% 
  select(new_id, everything())

#write_csv(ruv_k2_write, file = "DataProcessed/ruv_factor_data_k2_2022_04_20.csv")
