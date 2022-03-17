#load libraries
library(tidyverse)
library(readr)
library(RUVSeq)
library(here)
library(edgeR)
library(janitor)
library(DESeq2)

#read in raw gene counts
raw_gene_count <- read_tsv(file = here("DataRaw/gene_counts_choo.txt"), col_names = T)

#remove genes with all 0 counts


#Bioconductor suggested filter
filter <- raw_gene_count %>% 
  select(-Feature) %>% 
  apply(., 1, function(x) length(x[x>5])>=2)

filtered_gene_count <- raw_gene_count[filter,]
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

# #check Choo's original filter (0 count for more than 75% of samples)
# filtered_gene_count <- filtered_gene_count %>% 
#   rowwise() %>% 
#   mutate(sum_to_0_75 = sum(c_across(-Feature) == 0) >= 49*0.75)
# 
# #check how many genes that removes
# filtered_gene_count %>% 
#   group_by(sum_to_0_75) %>% 
#   summarise(N = n())
# 
# #filter them out
# filtered_gene_count <- filtered_gene_count %>% 
#   filter(sum_to_0_75 == F) %>% 
#   select(-sum_to_0_75)
# 
# #Check Choo's range filter
# filtered_gene_count <- filtered_gene_count %>% 
#   rowwise() %>% 
#   mutate(range_above_100 = (max(c_across(-Feature) - min(c_across(-Feature)))) >= 100)
# 
# #Check how many genes that removes
# filtered_gene_count %>% 
#   group_by(range_above_100) %>% 
#   summarise(N = n())
# 
# #Remove the genes
# filtered_gene_count <- filtered_gene_count %>% 
#   filter(range_above_100 == T) %>% 
#   select(-range_above_100)

##Prepare for RUV

#Load metadata
id_relate <- read_tsv(file = here("DataRaw/20201216_coreID_to_PID.txt"), col_names = T) %>% clean_names()
metadata_unjoined <- read_csv(file = here("DataProcessed/metadata_cleaning/table1_clean_data_2022_03_02.csv"))

#Join metadata
metadata_joined <- id_relate %>% 
  mutate(new_id = str_pad(new_id, 2, "left", "0") %>% paste0("Sample", .)) %>% 
  left_join(metadata_unjoined, by = "sid") %>% 
  filter(new_id %in% names(filtered_gene_count))

#Set up factors properly
metadata_joined$male_lab <- factor(metadata_joined$male_lab, 
                                   levels = c("Not Male", "Male"))
metadata_joined$vape_6mo_lab <- factor(metadata_joined$vape_6mo_lab, 
                                       levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"))
metadata_joined$latino_lab <- factor(metadata_joined$latino_lab, 
                                     levels = c("Non-LatinX", "LatinX"))

#Filter out the subject with missing vape status
metadata_joined <- metadata_joined %>% 
  filter(new_id != 'Sample23')

#make easier to reference
vape_status <- metadata_joined$vape_6mo_lab

male <- metadata_joined$male_lab

latinx <- metadata_joined$latino_lab

#convert to SEqExpressionSet object
ruv_prep <- newSeqExpressionSet(as.matrix(filtered_gene_count), 
                                phenoData = data.frame(vape_status, 
                                                       male, 
                                                       latinx, 
                                                       row.names = metadata_joined$new_id))


ruv_prep


#Look at raw sample
library(RColorBrewer)
library(ggpubr)
library(patchwork)
colors <- brewer.pal(3,"Set2")
par(mfrow = c(1,1))
plotRLE(ruv_prep, outline=FALSE, ylim=c(-4, 4), col = colors[metadata_joined$vape_6mo_lab])
plotPCA(ruv_prep, col = colors[metadata_joined$vape_6mo_lab])

# #considering upper-quartile normalization
# ruv_prep <- betweenLaneNormalization(ruv_prep, which = "upper")
# plotRLE(ruv_prep, outline=FALSE, ylim=c(-4, 4), col = colors[metadata_joined$vape_6mo_lab])

#Get the residuals matrix
#Design Matrix
design <- model.matrix(~ vape_status + male + latinx, data = pData(ruv_prep))
#get library size for each sample
ruv_counts <- DGEList(counts = counts(ruv_prep), group = vape_status)
#get commen negative binomial dispersion factor
ruv_counts <- estimateGLMCommonDisp(ruv_counts, design = design)


#fit log-linear negative binomial generalized model to the read counts for each gene
first_pass <- glmFit(ruv_counts, design = design)
#get the residuals from that model
first_pass_residuls <- residuals(first_pass, type="deviance")

## Run RUVr
#K = 1
ruv_k1 <- RUVr(ruv_prep, cIdx = genes, k = 1, residuals = first_pass_residuls)
pData(ruv_k1)

plotRLE(ruv_k1, outline=FALSE, ylim=c(-4, 4), col = colors[metadata_joined$vape_6mo_lab])

#K = 2
ruv_k2 <- RUVr(ruv_prep, cIdx = genes, k = 2, residuals = first_pass_residuls)

plotRLE(ruv_k2, outline=FALSE, ylim=c(-4, 4), col = colors[metadata_joined$vape_6mo_lab])

#K = 3
ruv_k3 <- RUVr(ruv_prep, cIdx = genes, k = 3, residuals = first_pass_residuls)

plotRLE(ruv_k3, outline=FALSE, ylim=c(-4, 4), col = colors[metadata_joined$vape_6mo_lab])

#K = 4
ruv_k4 <- RUVr(ruv_prep, cIdx = genes, k = 4, residuals = first_pass_residuls)

plotRLE(ruv_k4, outline=FALSE, ylim=c(-4, 4), col = colors[metadata_joined$vape_6mo_lab])

#K = 5
ruv_k5 <- RUVr(ruv_prep, cIdx = genes, k = 5, residuals = first_pass_residuls)

x <- plotRLE(ruv_k5, outline=FALSE, ylim=c(-4, 4), col = colors[metadata_joined$vape_6mo_lab])



