#Testing for correlation between center and vape status
#load libraries
library(tidyverse)
library(readr)
library(here)
library(vcd)
library(proxyC)
library(janitor)
library(rstatix)

###############################################################################################
##GET P-Value Histograms
#read in raw gene counts
raw_gene_count <- read_tsv(file = here("DataRaw/RNA_Seq/gene_counts_choo.txt"), col_names = T)

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
metadata_unjoined <- read_csv(file = here("DataProcessed/clinical_metadata/table1_clean_data_2022_08_22.csv"))


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
tab1_dat <- read_csv(file = here("DataProcessed/clinical_metadata/table1_clean_data_2022_08_22.csv"))

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

#Get p-values for correlation between each variable
for (i in variables) {
  for (j in variables) {
    Variable_1 <- base::append(Variable_1, i)
    Variable_2 <- base::append(Variable_2, j)
    #Create a contingency table
    cont_tab <- base::table(vars_of_interest_mat[,i],vars_of_interest_mat[,j])
    #if there are small counts in the matrix, run Fisher's. Else run pearson's chisq
    pval <- as.numeric(stats::fisher.test(cont_tab)[1])
    cor_test_mat[i,j] <- pval
    }
}


#Get upper-triangular of matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

cor_upper_tri <- get_upper_tri(cor_test_mat)

library(reshape2)
melted_cormat <- melt(cor_upper_tri, na.rm = T)

melted_cormat <- melted_cormat %>% 
  filter(Var1 != Var2)


correlation_matrix_plot <- melted_cormat %>% 
  ggplot(aes(x = Var2,y = Var1, fill = value)) + 
  geom_tile(color = "white") +
  scale_x_discrete(labels = c("Rectruitment Center", "LatinX", "Sex")) + 
  scale_y_discrete(labels = c("Vape Status", "Rectruitment Center", "LatinX")) +
  geom_text(aes(label = if_else(value != 1, prettyNum(value, preserve.width = "common", digits = 3), "> 0.99")), color = "white") +
  labs(fill = "P-Value")

correlation_matrix_plot

#complete t-tests to check for correlation between vape status 

age_t_test_tab <- tab1_dat %>% 
  select(sid,vape_6mo_lab, age) %>% 
  filter(sid %in% metadata_joined$sid) %>% 
  select(-sid) %>% 
  mutate(vape_status = if_else(vape_6mo_lab == "Vaped in Last 6 Months", 1, 0))

age_t_test <- t.test(age_t_test_tab$age~age_t_test_tab$vape_6mo_lab)


age_report <- tibble("Test Variable" = "Age",
                     "Group 1" = "Did Not Vape in Last 6 Months",
                     "Group 2" = "Vaped in last 6 months",
                     "Mean Group 1" = age_t_test$estimate[1],
                     "Mean Group 2" = age_t_test$estimate[2],
                     "t" = age_t_test$statistic,
                     "df" = age_t_test$parameter,
                     "p-value" = age_t_test$p.value
                     )
age_report_tab <- kable(age_report, digits = 1, booktabs = T) %>% 
  kable_styling(latex_options = "striped")
