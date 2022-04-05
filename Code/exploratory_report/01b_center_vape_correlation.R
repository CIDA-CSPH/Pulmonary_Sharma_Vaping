#Testing for correlation between center and vape status
#load libraries
library(tidyverse)
library(readr)
library(RUVSeq)
library(here)
library(edgeR)
library(janitor)
library(DESeq2)
library(dendextend)
library(WGCNA)

#read in raw gene counts
raw_gene_count <- read_tsv(file = here("DataRaw/gene_counts_choo.txt"), col_names = T)


#remove all genes with 0 counts in 75% or more of samples
filter_0_.75_count <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) sum(x == 0) >= 49*.75)


zero_count_.75_filtered <- raw_gene_count[!filter_0_.75_count,]

#remove all genes with range >= 100
filter_range_100 <- zero_count_.75_filtered %>% 
  select(-Feature) %>% 
  apply(.,1, function(x) (max(x) - min(x)) >= 100)

filtered_gene_count <- zero_count_.75_filtered[!filter_range_100,]

#convert gene count matrix
filtered_gene_count <- as.data.frame(filtered_gene_count)
rownames(filtered_gene_count) <- filtered_gene_count[,1]
filtered_gene_count <- filtered_gene_count[,-1]

#create a list of genes
genes <- rownames(filtered_gene_count)[grep("^ENS", rownames(filtered_gene_count))]

filtered_gene_count$Feature[1] == genes[1]

#Filter out Sample23(no vape status)
filtered_gene_count <- filtered_gene_count %>% 
  select(-Sample23)

##Prepare for RUV

#Load metadata
id_relate <- read_tsv(file = here("DataRaw/20201216_coreID_to_PID.txt"), col_names = T) %>% clean_names()
metadata_unjoined <- read_csv(file = here("DataProcessed/metadata_cleaning/table1_clean_data_2022_03_30.csv"))

#Join metadata
metadata_joined <- id_relate %>% 
  mutate(new_id = str_pad(new_id, 2, "left", "0") %>% paste0("Sample", .)) %>% 
  left_join(metadata_unjoined, by = "sid") %>% 
  filter(new_id %in% names(filtered_gene_count))

#Set up factors properly
metadata_joined$sex_lab <- factor(metadata_joined$sex_lab, 
                                  levels = c("Female", "Male"))
metadata_joined$vape_6mo_lab <- factor(metadata_joined$vape_6mo_lab, 
                                       levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"))
metadata_joined$latino_lab <- factor(metadata_joined$latino_lab, 
                                     levels = c("Non-LatinX", "LatinX"))

#Filter out the subject with missing vape status
metadata_joined <- metadata_joined %>% 
  filter(new_id != 'Sample23')

#make easier to reference
vape_status <- metadata_joined$vape_6mo_lab

male <- metadata_joined$sex_lab

latinx <- metadata_joined$latino_lab

#it appears that sample 12 may be outlier (see teen_vape_exploratory_report)

#Remove sample 12 as outlier
metadata_joined_no12 <- metadata_joined %>% 
  filter(new_id != 'Sample12')

filtered_gene_count_no12 <- filtered_gene_count %>% 
  select(-Sample12)


#make easier to reference
vape_status <- metadata_joined_no12$vape_6mo_lab

male <- metadata_joined_no12$sex_lab

latinx <- metadata_joined_no12$latino_lab

center <- metadata_joined_no12$recruitment_center

#check for correlation of vape_status and center

vape_center_cor_tab <- base::table(vape_status, center)

stats::fisher.test(vape_center_cor_tab)

#convert to SEqExpressionSet object
seq_mod <- newSeqExpressionSet(as.matrix(filtered_gene_count_no12), 
                                 phenoData = data.frame(vape_status, 
                                                        male, 
                                                        latinx,
                                                        center,
                                                        row.names = metadata_joined_no12$new_id))



##Get the residuals matrix
#Design Matrices
design_no_center <- model.matrix(~ vape_status + male + latinx, data = pData(seq_mod))
design_center <- model.matrix(~ vape_status + center + male + latinx, data = pData(seq_mod))
#get library size for each sample
counts <- DGEList(counts = counts(seq_mod), group = vape_status)
#get commen negative binomial dispersion factor
counts_no_center <- estimateGLMCommonDisp(counts, design = design_no_center)
counts_center <- estimateGLMCommonDisp(counts, design = design_center)
#fit model
fit_no_center <- glmQLFit(counts_no_center, design = design_no_center)
fit_center <- glmQLFit(counts_center, design = design_center)

#Test for DE Genes for model with and without center
test_no_center <- glmQLFTest(fit_no_center)
test_center <- glmQLFTest(fit_center)

test_pvals <- tibble(no_center = test_no_center$table$PValue,
                     center = test_center$table$PValue)
pvals_plot <- pivot_longer(test_pvals, cols = everything(), names_to = "center_included")

par(mfrow = c(1,2))
hist(test_pvals$no_center, xlab = "p-values", main = "No Center in Model", ylim = c(0,1000))
hist(test_pvals$center, xlab = "p-values", main = "Center in Model", ylim = c(0,1000))

pvals_plot %>% 
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~ center_included)
