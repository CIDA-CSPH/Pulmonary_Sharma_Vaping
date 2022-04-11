#Testing for correlation between center and vape status
#load libraries
library(tidyverse)
library(readr)
library(here)
library(vcd)
library(proxyC)

###############################################################################################
##GET P-Value Histograms
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

###############################################################################################
#Test for Correlation between all variables of Interest########################################
#Load variables 
tab1_dat <- read_csv(file = here("DataProcessed/metadata_cleaning/table1_clean_data_2022_03_30.csv"))

#Subset to only variables of interest
vars_of_interest <- tab1_dat %>% 
  select(sid,vape_6mo_lab, recruitment_center, latino_lab, sex_lab) %>% 
  filter(sid %in% metadata_joined$sid) %>% 
  select(-sid)

vars_of_interest_mat <- as.matrix(vars_of_interest)

variables <- base::colnames(vars_of_interest)

#create a results matrix
cor_test_mat <- matrix(nrow = ncol(vars_of_interest), ncol = ncol(vars_of_interest), 
                       dimnames = list(variables,variables))

Variable_1 <- vector()
Variable_2 <- vector()
test_type <- vector()

#Get p-values for correlation between each variable
for (i in variables) {
  for (j in variables) {
    Variable_1 <- base::append(Variable_1, i)
    Variable_2 <- base::append(Variable_2, j)
    #Create a contingency table
    cont_tab <- base::table(vars_of_interest_mat[,i],vars_of_interest_mat[,j])
    #if there are small counts in the matrix, run Fisher's. Else run pearson's chisq
    if (sum(apply(cont_tab, 2, function(x) sum(x < 5))) > 0) {
      pval <- as.numeric(stats::fisher.test(cont_tab)[1])
      cor_test_mat[i,j] <- pval
      test_type <- base::append(test_type, "Fisher Exact")
    }else{
      pval <- as.numeric(stats::chisq.test(cont_tab)[3])
      cor_test_mat[i,j] <- pval
      test_type <- base::append(test_type, "Pearson Chisq")
    }
    
  }
}
#Get test_type for each test
test_type_tracked <- as_tibble(cbind(Variable_1, Variable_2, test_type)) %>% 
  mutate(test_type_match = paste0(Variable_1, ", ", Variable_2))


#Get upper-triangular of matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

cor_upper_tri <- get_upper_tri(cor_test_mat)

library(reshape2)
melted_cormat <- melt(cor_upper_tri, na.rm = T)

#match up test_type for each test
melted_cormat <- melted_cormat %>% 
  dplyr::mutate(test_type_match = paste0(Var1, ", ", Var2))

test_type_matched <- test_type_tracked %>% 
  filter(test_type_match %in% melted_cormat$test_type_match)

melted_cormat <- left_join(melted_cormat, test_type_matched, by = "test_type_match")

#  
melted_cormat <- melted_cormat %>% 
  select(-c(Variable_1, Variable_2, test_type_match))

correlation_matrix_plot <- melted_cormat %>% 
  ggplot(aes(x = Var2,y = Var1, fill = if_else(value < 0.05,"<0.05",">0.05"))) + 
  geom_tile(color = "white") +
  scale_x_discrete(labels = c("Vape Status", "Rectruitment Center", "LatinX", "Sex")) + 
  scale_y_discrete(labels = c("Vape Status", "Rectruitment Center", "LatinX", "Sex")) +
  geom_text(aes(label = test_type), color = "white") +
  labs(fill = "Significance")
