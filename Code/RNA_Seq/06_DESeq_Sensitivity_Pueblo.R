######################### Load Libraries ###############################
library(ggpubr)
library(tidyverse)
library(here)
library(DESeq2)
library(RColorBrewer)
library(kableExtra)
library(ggbreak)

######################### Read in files so that they are unaltered and join ###############################
filtered_gene_count <- as.data.frame(read_csv(here("DataProcessed/rna_seq/differential_expression/full_analysis/filtered_gene_count_2022_10_13.csv")))
metadata_joined <- as.data.frame(read_csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv")))
ruv_factor_dat <- read_csv(here("DataProcessed/rna_seq/ruv/ruv_factor_data_k2_2022_10_13.csv"))
ruv_norm_counts <- as.data.frame(read_csv(here("DataProcessed/rna_seq/ruv/RUV_k2_norm_counts_2022_10_13.csv")))
gene_annotations <- read_tsv(here("DataRaw/RNA_Seq/gencode_annotations_choo.txt"))
de_full_res <- read_csv(here("DataProcessed/rna_seq/differential_expression/full_analysis/full_vape_res_2022_10_13.csv"))

#fix rownames of filtered gene count and ruv_counts
rownames(filtered_gene_count) <- filtered_gene_count$Feature
filtered_gene_count <- filtered_gene_count %>% 
  dplyr::select(-Feature)

#join metadata with ruv factor data
metadata_joined <- left_join(metadata_joined,ruv_factor_dat, by = "rna_id", copy = T) %>% 
  select(-c(vape_status, age.y, male)) 
metadata_joined <- metadata_joined %>%
  mutate(age = age.x, ruv_k1 = W_1, ruv_k2 = W_2) %>% 
  select(-c(age.x, W_1, W_2)) %>% 
  drop_na(ruv_k1)

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
rownames(metadata_joined) <- metadata_joined$rna_id


#check that row and column names are equal
all(rownames(metadata_joined) == colnames(filtered_gene_count))

######################### Filter for Pueblo ###############################
#metadata
metadata_joined_pblo <- metadata_joined %>%
  filter(recruitment_center == 'Pueblo')

 #count data
filtered_gene_count_pblo <- filtered_gene_count[,rownames(metadata_joined_pblo)]

#ruv normcounts
ruv_norm_counts_pblo <- ruv_norm_counts [,rownames(metadata_joined_pblo)]

######################### Design Matrices ###############################
full_design <- ~ sex_lab + age + ruv_k1 + ruv_k2 + vape_6mo_lab

vape_only_reduced <- ~ sex_lab + age + ruv_k1 + ruv_k2

######################### Required Functions ###############################
# format numbers
format_num <- function(number, digits = 0) {
  formatC(number, digits = digits, format = "f", big.mark = ",")
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
  deseq_res <- DESeq2::results(deseq_run, tidy = T, alpha = 0.05, name = "vape_6mo_lab_vaped_vs_not_vaped") %>% 
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

#tidy results for plotting (DESEQ-Normalized)
de_res_tidy <- function(run_de, de_res, col_data) {
  #Get top 4 genes of interest
  genes_of_interest <- de_res$row[1:4]
  
  #Normalize Counts
  run_de <- estimateSizeFactors(run_de)
  
  #Join and tidy
  norm_counts<- t(log10((counts(run_de[genes_of_interest, ], 
                                normalized=TRUE)+.5))) %>%
    merge(col_data, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(genes_of_interest) + 1):ncol(.))
  #add symbol name
  gene_symbols <- gene_annotations[gene_annotations$gene %in% genes_of_interest,]
  
  norm_counts <- left_join(norm_counts, gene_symbols, by = "gene") 
  
  return(norm_counts)
}

#tidy results for plotting (RUV-Normalized)
ruv_res_tidy <- function(ruv_count_dat, de_res, col_data) {
  #Get top 4 genes of interest
  genes_of_interest <- de_res$row[1:4]
  
  #Join and tidy
  norm_counts<- t(log10((ruv_count_dat[genes_of_interest,] +.5))) %>%
    merge(col_data, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(genes_of_interest) + 1):ncol(.))
  #add symbol name
  gene_symbols <- gene_annotations[gene_annotations$gene %in% genes_of_interest,]
  
  norm_counts <- left_join(norm_counts, gene_symbols, by = "gene") 
  
  return(norm_counts)
}

#Results Boxplots (Vape and Center, Vape Only, Center Only)
tcount_boxplot <- function(tcounts, mod) {
  # #Stratified by Vape and Center
  # vape_center_box <- tcounts %>% 
  #   ggplot(aes(vape_6mo_lab, expression, group = interaction(vape_6mo_lab, recruitment_center))) + 
  #   geom_boxplot(aes(fill = recruitment_center)) +
  #   geom_text(label = if_else(tcounts$Row.names == "Sample12", tcounts$Row.names, ""), col = "Red") +
  #   facet_wrap(~symbol, scales="free_y") + 
  #   labs(x="", 
  #        y="Expression (log normalized counts)", 
  #        fill="Center", 
  #        title= paste0("Top Genes (", mod, ")")) +
  #   scale_fill_manual(values = center_colors) + 
  #   theme(legend.position = "bottom")
  
  #Stratified by vape status only
  vape_box <- tcounts %>% 
    ggplot(aes(vape_6mo_lab, expression, fill= vape_6mo_lab)) + 
    geom_boxplot() + 
    # geom_text(label = if_else(tcounts$Row.names == "Sample12", tcounts$Row.names, ""), col = "Red") +
    labs(x="Vape Status (6 mo)", 
         y="Expression (log normalized counts)", 
         fill = "Vape Status (6 mo)", 
         title= paste0("Top Genes (", mod, ")")) +
    scale_fill_manual(values = vape_colors) +
    facet_wrap(~symbol, scales="free_y") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom")
  
  # #Stratified by center only
  # center_box <- tcounts %>% 
  #   ggplot(aes(recruitment_center, expression, fill=recruitment_center)) + 
  #   geom_boxplot() +
  #   geom_text(label = if_else(tcounts$Row.names == "Sample12", tcounts$Row.names, ""), 
  #             col = "Red") +
  #   facet_wrap(~symbol, scales="free_y") + 
  #   labs(x="Recruitment Center", 
  #        y="Expression (log normalized counts)", 
  #        fill="Center", 
  #        title= paste0("Top Genes (", mod, ")")) +
  #   scale_fill_manual(values = center_colors) + 
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         legend.position = "bottom") 
  
  return(vape_box)
  
}

#Counts Significant Genes from an object created using de_res_tidy function
sig_gene_count <- function(de_res) {
  sum(de_res$padj < 0.05)
}

######################### Run DESeq ###############################
#Vape Only
vape_de_pblo <- run_deseq_lrt(count_data = filtered_gene_count_pblo, 
                              col_data = metadata_joined_pblo,
                              full_mod = full_design,
                              reduced_mod = vape_only_reduced)

######################### Get Formatted Results ###############################

vape_res_pblo <- format_results(vape_de_pblo, gene_annotations)

#Write out results
write_csv(vape_res_pblo, "DataProcessed/rna_seq/differential_expression/pueblo_sensitivity/vape_res_pblo_2022_10_21.csv")

######################### Filter for log-fold change (check) ###############################

check_lf_cutoff <- function(de_res) {
  cutoff_seq <- seq(0,6,0.5)
  cutoff_res <- NULL
  for (cutoff in cutoff_seq) {
    num_genes <- nrow(de_res %>%
                        filter(padj < 0.05) %>% 
                        filter(abs(log2FoldChange) > cutoff))
    cutoff_res <- append(cutoff_res, num_genes)
  }
  
  result <- tibble(Cutoff = cutoff_seq,
                   Genes = cutoff_res)
  return(result)
}

pblo_res_cutoff <- check_lf_cutoff(vape_res_pblo)

full_res_cutoff <- check_lf_cutoff(de_full_res)

full_res_cutoff %>% 
  ggplot(aes(x = Cutoff, y = Genes)) +
  geom_point() +
  geom_line() + geom_hline(yintercept = 1000, col = 'red')

######################### Implement log2 Fold Change cutoff ###############################

#pueblo sensitivity
vape_res_pblo <- vape_res_pblo %>% 
  filter(abs(log2FoldChange) > 1)

#full analysis
de_full_res <- de_full_res %>% 
  filter(abs(log2FoldChange) > 1)

######################### P-Value Histograms ###############################
#Vape Only
vape_only_phist_pblo <- p_hist(vape_res_pblo) +
  labs(title = "Vape Only")

vape_only_phist_pblo
######################### Top Genes Tidy DE-Norm Counts ###############################
#fix metadata_joined vape status
metadata_joined_pblo$vape_6mo_lab <- if_else(metadata_joined_pblo$vape_6mo_lab == "vaped", "Vaped", "Not Vaped")

#Vape only
vape_de_tidy_pblo  <- de_res_tidy(de_res = vape_res_pblo, 
                                  run_de = vape_de_pblo,
                                  col_data = metadata_joined_pblo)

######################### Top Genes Tidy RUV-Norm Counts ###############################

vape_ruv_tidy_pblo <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts_pblo,
                                   de_res = vape_res_pblo,
                                   col_data = metadata_joined_pblo)



# ######################### Top Gene Boxplots DE Counts ###############################
# #Set color pallettes to match
# vape_colors <- brewer.pal(3,"Set2")
# 
# #Vape only
# vape_plots_pblo <- tcount_boxplot(tcounts = vape_de_tidy_pblo, mod = "Vape Only")
# 
# 
# vape_plots_pblo
# 
# 
# ######################### Top Gene Boxplots RUV Counts ###############################
# 
# #Vape only
# vape_plots_ruv_pblo <- tcount_boxplot(tcounts = vape_ruv_tidy_pblo, mod = "Vape Only")
# 
# 
# vape_plots_ruv_pblo

