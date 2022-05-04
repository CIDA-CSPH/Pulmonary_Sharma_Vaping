#load libraries
library(tidyverse)
library(readr)
library(janitor)
library(RUVSeq)
library(here)
library(DESeq2)
library(Hmisc)
library(reshape2)
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

genes_filtered <- rownames(filtered_gene_count)

##Prepare for DESeq2
#Set up factors properly
metadata_joined$sex_lab <- factor(metadata_joined$sex_lab, 
                                  levels = c("Female", "Male"))
metadata_joined$vape_6mo_lab <- factor(metadata_joined$vape_6mo_lab, 
                                       levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"),
                                       labels = c("Did_Not_Vape_in_Last_6_Months", "Vaped_in_Last_6_Months"))
metadata_joined$latino_lab <- factor(metadata_joined$latino_lab, 
                                     levels = c("Non-LatinX", "LatinX"),
                                     labels = c("Non_LatinX", "LatinX"))


#Assign sample names to rownames
rownames(metadata_joined) <- metadata_joined$new_id
metadata_joined <- metadata_joined %>% 
  select(-new_id)

#check that row and column names are equal
all(rownames(metadata_joined) == colnames(filtered_gene_count))

#make design matrix
design_ruv <- ~vape_6mo_lab + sex_lab + latino_lab

ruv_prep <- DESeqDataSetFromMatrix(countData = filtered_gene_count,
                                   colData = metadata_joined,
                                   design = design_ruv) 
#Run DESeq First Pass

DESeq_first_pass_vst <- DESeq(ruv_prep) %>% vst() %>% assay()

DESeq_Resid <-
  DESeq_first_pass_vst %>%
  apply(., 1, function(y) {
    lm(y ~ vape_6mo_lab +
         latino_lab + sex_lab,
       data = metadata_joined) %>% resid()
  })

DESeq_Resid <- DESeq_Resid %>% t()

#write out needed objects
#write_csv(as.data.frame(DESeq_Resid), file = here("DataProcessed/first_pass_residuals_DESeq2.csv"))
#write_csv(metadata_joined, file = "DataProcessed/metadata_joined_04_20_2022.csv")
# write_csv(filtered_gene_count, file = "DataProcessed/filtered_gene_count_04_20_2022.csv")
# write(genes_filtered, file = "DataProcessed/filtered_gene_list_04_20_2022.txt")
