#load libraries
library(tidyverse)
library(readr)
library(RUVSeq)
library(here)
library(janitor)
library(edgeR)
library(dendextend)
library(WGCNA)
library(Hmisc)

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
id_relate <- read_tsv(file = here("DataRaw/subject_ids/20201216_coreID_to_PID.txt"), col_names = T) %>% clean_names()
metadata_unjoined <- read_csv(file = here("DataProcessed/metadata/table1_clean_data_2022_08_22.csv"))


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
metadata_joined$latino_lab <- factor(metadata_joined$latino_lab, 
                                     levels = c("Non-LatinX", "LatinX"),
                                     labels = c("Non_LatinX", "LatinX"))



#check that row and column names are equal
#metadata_joined$new_id == colnames(filtered_gene_count)

##Prepare for RUV
#make easier to reference
vape_status <- metadata_joined$vape_6mo_lab

male <- metadata_joined$sex_lab

latinx <- metadata_joined$latino_lab

age <- metadata_joined$age

center <- metadata_joined$recruitment_center

#convert to SEqExpressionSet object WITHOUT Center
ruv_with_12 <- newSeqExpressionSet(as.matrix(filtered_gene_count), 
                                phenoData = data.frame(vape_status, 
                                                       male, 
                                                       age, 
                                                       row.names = metadata_joined$new_id))


#it appears that sample 12 may be outlier (see teen_vape_exploratory_report) will run RUV both with and without sample

#Remove sample 12 as outlier
metadata_joined_no12 <- metadata_joined %>% 
  filter(new_id != 'Sample12')

filtered_gene_count_no12 <- filtered_gene_count %>% 
  select(-Sample12)

#make easier to reference
vape_status_no12 <- metadata_joined_no12$vape_6mo_lab

male_no12 <- metadata_joined_no12$sex_lab

latinx_no12 <- metadata_joined_no12$latino_lab

age_no12 <- metadata_joined_no12$age

center_no12 <- metadata_joined_no12$recruitment_center

#convert to SEqExpressionSet object
ruv_no12 <- newSeqExpressionSet(as.matrix(filtered_gene_count_no12), 
                                phenoData = data.frame(vape_status_no12, 
                                                       male_no12, 
                                                       age_no12, 
                                                       row.names = metadata_joined_no12$new_id))


# #considering upper-quartile normalization
# ruv_prep <- betweenLaneNormalization(ruv_prep, which = "upper")
# plotRLE(ruv_prep, outline=FALSE, ylim=c(-4, 4), col = colors[metadata_joined$vape_6mo_lab])

##Get the residuals matrix
#Design Matrix
design <- model.matrix(~ vape_status + male + age, data = pData(ruv_with_12))
design_no12 <- model.matrix(~ vape_status_no12 + male_no12 + age_no12, data = pData(ruv_no12))
#get library size for each sample
ruv_counts <- DGEList(counts = counts(ruv_with_12), group = vape_status)
ruv_counts_no12 <- DGEList(counts = counts(ruv_no12), group = vape_status_no12)
#get commen negative binomial dispersion factor
ruv_counts <- estimateGLMCommonDisp(ruv_counts, design = design)
ruv_counts_no12 <- estimateGLMCommonDisp(ruv_counts_no12, design = design_no12)
#fit log-linear negative binomial generalized model to the read counts for each gene
first_pass <- glmFit(ruv_counts, design = design)
first_pass_no12 <- glmFit(ruv_counts_no12, design = design_no12)
#get the residuals from that model
first_pass_residuls <- residuals(first_pass, type="deviance")
first_pass_residuls_no12 <- residuals(first_pass_no12, type="deviance")

#output the firt pass residuals for comparison with DESeq2
#write_csv(as.data.frame(first_pass_residuls), file = here("DataProcessed/first_pass_residuals_edgeR.csv"))

#Elbow Method for finding the optimal number of clusters
set.seed(404)

#get scaled data
scaled_data <- as.matrix(scale(assayData(ruv_with_12)$counts))
scaled_data_no12 <- as.matrix(scale(assayData(ruv_no12)$counts))

# Compute and plot wss for k = 2 to k = 15.
k.max <- 5
wss <- sapply(1:k.max, 
              function(k){kmeans(scaled_data, k, nstart = 10, iter.max = 15)$tot.withinss})
wss_no12 <- sapply(1:k.max, 
                   function(k){kmeans(scaled_data_no12, k, nstart = 10, iter.max = 15)$tot.withinss})

elbow_tib <- tibble(k = seq(1,5,1), 
                    wss = wss,
                    wss_no12 = wss_no12)

elbow_tib_long <- pivot_longer(elbow_tib, cols = c(wss, wss_no12), names_to = "Sample") %>% 
  mutate(Sample = if_else(Sample == "wss", "With Sample 12", "Without Sample 12"))

elbow_plot <- elbow_tib_long %>% 
  ggplot(aes(x = k, y = value, col = Sample)) +
  geom_point() +
  geom_line() +
  geom_point(aes(x = ))+
  labs(x = "K", y = "Total Within-Clusters Sum of Squares", title = "Elbow Plot", col = "") + 
  scale_x_continuous(breaks = seq(1,10,1))
  

## Run RUVr
#K = 1
ruv_k1 <- RUVr(ruv_with_12, cIdx = genes, k = 1, residuals = first_pass_residuls)
ruv_k1_no12 <- RUVr(ruv_no12, cIdx = genes, k = 1, residuals = first_pass_residuls_no12)

#K = 2
ruv_k2 <- RUVr(ruv_with_12, cIdx = genes, k = 2, residuals = first_pass_residuls)
ruv_k2_no12 <- RUVr(ruv_no12, cIdx = genes, k = 2, residuals = first_pass_residuls_no12)

#K = 3
ruv_k3 <- RUVr(ruv_with_12, cIdx = genes, k = 3, residuals = first_pass_residuls)
ruv_k3_no12 <- RUVr(ruv_no12, cIdx = genes, k = 3, residuals = first_pass_residuls_no12)


