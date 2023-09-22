######################### Load Libraries ###############################
library(ggpubr)
library(tidyverse)
library(here)
library(DESeq2)
library(RColorBrewer)
library(kableExtra)

######################### Read in files so that they are unaltered and join ###############################
filtered_gene_count <- as.data.frame(read_csv(here("DataProcessed/rna_seq/differential_expression/full_analysis/filtered_gene_count_2022_05_04.csv")))
metadata_joined <- as.data.frame(read_csv(here("DataProcessed/clinical_metadata/metadata_joined_rnaseq_04_20_2022.csv")))
ruv_factor_dat <- read_csv(here("DataProcessed/rna_seq/ruv/ruv_factor_data_k2_2022_04_20.csv"))
ruv_norm_counts <- as.data.frame(read_csv(here("DataProcessed/rna_seq/ruv/RUV_k2_norm_counts_2022_05_06.csv")))
gene_annotations <- read_tsv(here("DataRaw/RNA_Seq/gencode_annotations_choo.txt"))

#fix rownames of filtered gene count and ruv_counts
rownames(filtered_gene_count) <- filtered_gene_count$Feature
filtered_gene_count <- filtered_gene_count %>% 
  dplyr::select(-Feature)

#join metadata with ruv factor data
metadata_joined <- left_join(metadata_joined,ruv_factor_dat, by = "new_id", copy = T) %>% 
  select(-c(vape_status, age.y, male)) 
metadata_joined <- metadata_joined %>%
  mutate(age = age.x, ruv_k1 = W_1, ruv_k2 = W_2) %>% 
  select(-c(age.x, W_1, W_2))

#fix gene annotations
gene_annotations <- gene_annotations %>% 
  dplyr::rename(row = ENSG)

#fix ruv colnames
rownames(ruv_norm_counts) <- ruv_norm_counts$gene
ruv_norm_counts <- ruv_norm_counts %>% 
  select(-gene)

######################### Prepare Metadata ###############################
#Set factor levels
metadata_joined$vape_6mo_lab <- factor(metadata_joined$vape_6mo_lab, levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"), 
                                       labels = c("not_vaped", "vaped"))

metadata_joined$sex_lab <- factor(metadata_joined$sex_lab, levels = c("Female", "Male"))

metadata_joined$recruitment_center <- factor(metadata_joined$recruitment_center, levels = c("Aurora", "CommCity/Denver", "Pueblo"), 
                                             labels = c("Aurora", "CommCity_Denver", "Pueblo"))
#Scaling Age
metadata_joined$age <- scale(metadata_joined$age)

#Make Sample ID into row names
rownames(metadata_joined) <- metadata_joined$new_id

metadata_joined <- metadata_joined %>% 
  select(-new_id)

#check that row and column names are equal
all(rownames(metadata_joined) == colnames(filtered_gene_count))

######################### Design Matrices ###############################
full_design <- ~vape_6mo_lab + sex_lab + age + recruitment_center + ruv_k1 + ruv_k2 

reduced <- ~ sex_lab + age + recruitment_center + ruv_k1 + ruv_k2

######################### Required Functions ###############################
# format numbers
format_num <- function(number, digits = 0) {
  formatC(x, digits = digits, format = "f", big.mark = ",")
}

#Run DESeq
run_deseq_lrt <- function(count_data, col_data, full_mod, reduced_mod) {
  
  #Create the DESeq Object
  deseq_object <- DESeqDataSetFromMatrix(countData = count_data,
                                         colData = col_data,
                                         design = full_mod)
  #Run DESeq
  deseq_run <- DESeq(deseq_object, test = "LRT", reduced = reduced_mod)
  
  #Return the results
  return(deseq_run)
}

#DESeq Results table
format_results <- function(deseq_run, annotation) {
  #Pull out resultsand sort by padj and pval
  deseq_res <- DESeq2::results(deseq_run, tidy = T, alpha = 0.05) %>% 
    arrange(padj, pvalue) %>% 
    as_tibble() %>% 
    left_join(.,annotation, by = 'row') %>% 
    select(row, symbol, gene_type, everything()) %>% 
    dplyr::rename(ensg = row,
                  gene_name = symbol)
  return(deseq_res)
}

#p-value histograms
p_hist <- function(de_results) {
  temp_hist <- de_results %>% 
    ggplot(aes(x = pvalue))+
    geom_histogram()
  return(temp_hist)
}

#tidy results for plotting (RUV-Normalized)
ruv_res_tidy <- function(ruv_count_dat, de_res, col_data, annotation) {
  #Get top 4 genes of interest
  genes_of_interest <- de_res$ensg[1:4]
  
  #Join and tidy
  norm_counts<- t(log10((ruv_count_dat[genes_of_interest,] +.5))) %>%
    merge(col_data, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(genes_of_interest) + 1):ncol(.))
  
  #join back with gene name
  norm_counts <- left_join(norm_counts, annotation, by = c("gene" = "row"))
  
  return(norm_counts)
}

#Results Boxplots (Vape and Center, Vape Only, Center Only)
tcount_boxplot <- function(tcounts, mod) {
  
  #Stratified by vape status only
  vape_box <- tcounts %>% 
    ggplot(aes(vape_6mo_lab, expression, fill= vape_6mo_lab)) + 
    geom_boxplot() + 
    geom_text(label = if_else(tcounts$Row.names == "Sample12", tcounts$Row.names, ""), col = "Red") +
    labs(x="Vape Status (6 mo)", 
         y="Expression (log normalized counts)", 
         fill = "Vape Status (6 mo)", 
         title= paste0("Top Genes (", mod, ")")) +
    scale_fill_manual(values = vape_colors) +
    facet_wrap(~symbol, scales="free_y") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom")
  
  return(vape_box)
  
}

#Counts Significant Genes from an object created using de_res_tidy function
sig_gene_count <- function(de_res) {
  sum(de_res$padj < 0.05)
}

######################### Run DESeq ###############################
#Vape Only
vape_de <- run_deseq_lrt(count_data = filtered_gene_count, 
                         col_data = metadata_joined,
                         full_mod = full_design,
                         reduced_mod = reduced)

######################### Get Formatted Results ###############################

vape_res <- format_results(vape_de, gene_annotations)

#Filter it for boxplots 
vape_res <- vape_res %>% 
  filter(abs(log2FoldChange) > 2)
#Write out results
# write_csv(vape_res, file = here("DataProcessed/de_full/full_vape_res_2022_05_01.csv"))

######################### P-Value Histograms ###############################

#Vape Only
vape_only_phist <- p_hist(vape_res) +
  labs(title = "Vape Only")

######################### Top Genes Tidy RUV-Norm Counts ###############################

vape_ruv_tidy <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts,
                              de_res = vape_res,
                              col_data = metadata_joined,
                              annotation = gene_annotations)

#write out results
# write_csv(vape_ruv_tidy, file = here("DataProcessed/de_full/top4_vape_ruv_normcounts_2022_06_14.csv"))


######################### Top Gene Boxplots RUV Counts ###############################
#Set color pallettes to match
vape_colors <- brewer.pal(3,"Set2")

#Vape only
vape_plots_ruv <- tcount_boxplot(tcounts = vape_ruv_tidy, mod = "Vape Only")


vape_plots_ruv

