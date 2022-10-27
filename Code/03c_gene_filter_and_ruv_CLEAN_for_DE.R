#load libraries
library(tidyverse)
library(readr)
library(RUVSeq)
library(here)
library(edgeR)
library(janitor)


#read in raw gene counts
raw_gene_count <- read.table(file = here("DataRaw/RNA_Seq/gene_counts_choo.txt"), row.names = 1, header = T) 

#remove all genes with 0 counts in 75% or more of samples
filter_0_.75_count <- raw_gene_count %>%
  apply(., 1, function(x) sum(x == 0) >= 49*.75)

#True if # of cells with 0's >= 75% of samples
#Keep if NOT true
zero_count_.75_filtered <- raw_gene_count[!filter_0_.75_count,]

#remove all genes with range >= 100
filter_range_100 <- zero_count_.75_filtered %>% 
  apply(.,1, function(x) (max(x) - min(x)) >= 100)

#True if range of read counts across all samples are >= 100
#Keep if TRUE. Want to remove genes with low variation across samples
filtered_gene_count <- zero_count_.75_filtered[filter_range_100,]

#create a list of genes
genes <- rownames(filtered_gene_count)

#Load metadata
metadata_joined <- as.data.frame(read_csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv")))

metadata_joined <- metadata_joined %>% 
  drop_na(rna_id, vape_6mo_lab) %>% 
  filter(rna_id != "Sample35")

#Make sure the count data matches
filtered_gene_count <- filtered_gene_count[,names(filtered_gene_count) %in% metadata_joined$rna_id]

filtered_gene_count_write <- filtered_gene_count %>% 
  rownames_to_column(var = "Feature")

#write_csv(filtered_gene_count_write, file = here("DataProcessed/filtered_gene_count_2022_10_13.csv"))

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
                                                          row.names = metadata_joined$rna_id))


#Read in residuals matrix (See "03_gene_filter_and_ruv_edgeR.R" for more information on calculation)
#Get the residuals matrix
#Design Matrix
design <- model.matrix(~ vape_status + male + age, data = pData(ruv_prep))

#get library size for each sample
ruv_counts <- DGEList(counts = counts(ruv_prep), group = vape_status)

#get commen negative binomial dispersion factor
ruv_counts <- estimateGLMCommonDisp(ruv_counts, design = design)

#fit log-linear negative binomial generalized model to the read counts for each gene
first_pass <- glmFit(ruv_counts, design = design)

#get the residuals from that model
first_pass_residuls <- residuals(first_pass, type="deviance")

rownames(first_pass_residuls) <- rownames(filtered_gene_count)
colnames(first_pass_residuls) <- colnames(filtered_gene_count)

## Run RUVr with k = 2
ruv_k2 <- RUVr(ruv_prep, cIdx = genes, k = 2, residuals = first_pass_residuls)

#write out needed results
ruv_k2_norm_counts <- as.data.frame(normCounts(ruv_k2)) %>% 
  rownames_to_column(var = "Feature")

#write_csv(ruv_k2_norm_counts, file = here("DataProcessed/rna_seq/ruv/RUV_k2_norm_counts_2022_10_13.csv"))

ruv_k2_write <- pData(ruv_k2) %>% 
  rownames_to_column(var = "rna_id")

#write_csv(ruv_k2_write, file = "DataProcessed/rna_seq/ruv/ruv_factor_data_k2_2022_10_13.csv")
