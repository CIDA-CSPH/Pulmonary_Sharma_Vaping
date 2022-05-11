######################### Load Libraries ###############################
library(ggpubr)
library(tidyverse)
library(here)
library(DESeq2)
library(pcaExplorer)
library(RColorBrewer)
library(kableExtra)
library(ggbreak)

######################### PLEASE RUN "04_DESeq2_DE_Model_fit.R" first to have tables WITH sample 12 for Comparison ###############################
source(here("Code/04_DESeq2_DE_Model_fit.R"))

######################### Read in files so that they are unaltered and join ###############################
filtered_gene_count <- as.data.frame(read_csv(here("DataProcessed/filtered_gene_count_2022_05_04.csv")))
genes <- scan(file = here("DataProcessed/filtered_gene_list_2022_04_20.txt"), character(), quote = "")
metadata_joined <- as.data.frame(read_csv(here("DataProcessed/metadata_joined_2022_04_20.csv")))
ruv_factor_dat <- read_csv(here("DataProcessed/ruv_factor_data_k2_2022_04_20.csv"))
ruv_norm_counts <- as.data.frame(read_csv(here("DataProcessed/RUV_k2_norm_counts_2022_05_06.csv")))
gene_annotations <- read_tsv(here("DataRaw/gencode_annotations.txt"))

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
  dplyr::mutate(gene = ENSG) %>% 
  dplyr::select(gene, symbol)

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

######################### Filter out Sample 12 ###############################
#metadata
metadata_joined_no12 <- metadata_joined %>%
  filter(sid != 102)
#count data
filtered_gene_count_no12 <- filtered_gene_count %>%
  select(-Sample12)
#ruv normcounts
ruv_norm_counts_no12 <- ruv_norm_counts %>% 
  select(-Sample12)

######################### Design Matrices ###############################
full_design <- ~vape_6mo_lab + sex_lab + age + recruitment_center + ruv_k1 + ruv_k2 

vape_only_reduced <- ~vape_6mo_lab + sex_lab + age + ruv_k1 + ruv_k2

center_only_reduced <- ~ sex_lab + age + recruitment_center + ruv_k1 + ruv_k2

reduced_design <- ~sex_lab + age + ruv_k1 + ruv_k2 
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
format_results <- function(deseq_run) {
  #Pull out resultsand sort by padj and pval
  deseq_res <- DESeq2::results(deseq_run, tidy = T, alpha = 0.05) %>% 
    arrange(padj, pvalue) %>% 
    tbl_df()
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
  #Stratified by Vape and Center
  vape_center_box <- tcounts %>% 
    ggplot(aes(vape_6mo_lab, expression, group = interaction(vape_6mo_lab, recruitment_center))) + 
    geom_boxplot(aes(fill = recruitment_center)) +
    geom_text(label = if_else(tcounts$Row.names == "Sample12", tcounts$Row.names, ""), col = "Red") +
    facet_wrap(~symbol, scales="free_y") + 
    labs(x="", 
         y="Expression (log normalized counts)", 
         fill="Center", 
         title= paste0("Top Genes (", mod, ")")) +
    scale_fill_manual(values = center_colors) + 
    theme(legend.position = "bottom")
  
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
  
  #Stratified by center only
  center_box <- tcounts %>% 
    ggplot(aes(recruitment_center, expression, fill=recruitment_center)) + 
    geom_boxplot() +
    geom_text(label = if_else(tcounts$Row.names == "Sample12", tcounts$Row.names, ""), 
              col = "Red") +
    facet_wrap(~symbol, scales="free_y") + 
    labs(x="Recruitment Center", 
         y="Expression (log normalized counts)", 
         fill="Center", 
         title= paste0("Top Genes (", mod, ")")) +
    scale_fill_manual(values = center_colors) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") 
  
  return(list(vape_center_box, vape_box, center_box))
  
}

#Counts Significant Genes from an object created using de_res_tidy function
sig_gene_count <- function(de_res) {
  sum(de_res$padj < 0.05)
}

######################### Run DESeq ###############################
# Vape and Center
vape_center_de_no12 <- run_deseq_lrt(count_data = filtered_gene_count_no12, 
                                col_data = metadata_joined_no12,
                                full_mod = full_design,
                                reduced_mod = reduced_design) 
#Vape Only
vape_de_no12 <- run_deseq_lrt(count_data = filtered_gene_count_no12, 
                         col_data = metadata_joined_no12,
                         full_mod = full_design,
                         reduced_mod = center_only_reduced)

#Center Only
center_de_no12 <- run_deseq_lrt(count_data = filtered_gene_count_no12, 
                           col_data = metadata_joined_no12,
                           full_mod = full_design,
                           reduced_mod = vape_only_reduced)

######################### Get Formatted Results ###############################
vape_center_res_no12 <- format_results(vape_center_de_no12)

vape_res_no12 <- format_results(vape_de_no12)

center_res_no12 <- format_results(center_de_no12)

######################### P-Value Histograms ###############################
#Vape and Center
vape_center_phist_no12 <- p_hist(vape_center_res_no12) +
  labs(title = "Vape and Center")
#Vape Only
vape_only_phist_no12 <- p_hist(vape_res_no12) +
  labs(title = "Vape Only")
#Center Only
center_only_phist_no12 <- p_hist(center_res_no12) + 
  labs(title = "Center Only")

######################### Top Genes Tidy DE-Norm Counts ###############################
#fix metadata_joined vape status
metadata_joined_no12$vape_6mo_lab <- if_else(metadata_joined$vape_6mo_lab == "vaped", "Vaped", "Not Vaped")

#Vape and Center
vape_center_de_tidy_no12 <- de_res_tidy(de_res = vape_center_res_no12, 
                                   run_de = vape_center_de_no12,
                                   col_data = metadata_joined_no12)
#Vape only
vape_de_tidy_no12  <- de_res_tidy(de_res = vape_res_no12, 
                            run_de = vape_de_no12,
                            col_data = metadata_joined_no12)
#Center Only
center_de_tidy_no12  <- de_res_tidy(de_res = center_res_no12, 
                              run_de = center_de_no12,
                              col_data = metadata_joined_no12)

######################### Top Genes Tidy RUV-Norm Counts ###############################
vape_center_ruv_tidy_no12 <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts_no12,
                                     de_res = vape_center_res_no12,
                                     col_data = metadata_joined_no12)
vape_ruv_tidy_no12 <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts_no12,
                                   de_res = vape_res_no12,
                                   col_data = metadata_joined_no12)

center_ruv_tidy_no12 <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts_no12,
                                     de_res = center_res_no12,
                                     col_data = metadata_joined_no12)

######################### Top Gene Boxplots DE Counts ###############################
#Set color pallettes to match
center_colors <- brewer.pal(10,"Paired")[4:6]
vape_colors <- brewer.pal(3,"Set2")

#Vape and Center
vape_center_plots_no12 <- tcount_boxplot(tcounts = vape_center_de_tidy_no12, mod = "Vape and Center")

#Vape only
vape_plots_no12 <- tcount_boxplot(tcounts = vape_de_tidy_no12, mod = "Vape Only")

#Center Only
center_plots_no12 <- tcount_boxplot(tcounts = center_de_tidy_no12, mod = "Center Only")

#Looking at plots
#Vape and Center
vape_center_plots_no12[[1]]
vape_plots_no12[[1]]
center_plots_no12[[1]]

#Vape only
vape_center_plots_no12[[2]]
vape_plots_no12[[2]]
center_plots_no12[[2]]

#Center Only
vape_center_plots_no12[[3]]
vape_plots_no12[[3]]
center_plots_no12[[3]]

######################### Top Gene Boxplots RUV Counts ###############################
#Vape and Center
vape_center_plots_ruv_no12 <- tcount_boxplot(tcounts = vape_center_ruv_tidy_no12, mod = "Vape and Center")

#Vape only
vape_plots_ruv_no12 <- tcount_boxplot(tcounts = vape_ruv_tidy_no12, mod = "Vape Only")

#Center Only
center_plots_ruv_no12 <- tcount_boxplot(tcounts = center_ruv_tidy_no12, mod = "Center Only")


#Vape and Center
vape_center_plots_ruv_no12[[1]]
vape_plots_ruv_no12[[1]]
center_plots_ruv_no12[[1]]

#Vape only
vape_center_plots_ruv_no12[[2]]
vape_plots_ruv_no12[[2]]
center_plots_ruv_no12[[2]]

#Center Only
vape_center_plots_ruv_no12[[3]]
vape_plots_ruv_no12[[3]]
center_plots_ruv_no12[[3]]

######################### Results Table FDR ###############################
sig_genes_tab <- tibble(Model = c("With Sample 12", "Sample 12 Removed"),
                        "FDR < 0.05" = c(format_num(sig_gene_count(de_res = vape_res)),
                                                           format_num(sig_gene_count(de_res = vape_res_no12)))) %>% 
                                           kbl(digits = 3)
sig_genes_tab



######################### Log-fold change ###############################
# add in symbols
# Subset the top 100 genes
top_100_12 <- vape_res[1:100,] %>% arrange(row)
top_100_no12 <- vape_res_no12[vape_res_no12$row %in% top_100_12$row,] %>% arrange(row)

top_100_vape_both <- left_join(top_100_12, top_100_no12, by = "row", suffix = c(".12",".no12"))

top_100_vape_both <- top_100_vape_both %>% 
  dplyr::mutate(fold_percent_change = ((log2FoldChange.12 - log2FoldChange.no12)/log2FoldChange.12)*100,
                fold_change = (log2FoldChange.12 - log2FoldChange.no12)) %>% 
  left_join(., gene_annotations, by = c("row" = "gene"))

p_change_plot <- top_100_vape_both %>% 
  ggplot(aes(x = fold_percent_change)) +
  geom_histogram(bins = 50, color="white") +
  labs(x = "% Change",
       y = "Count") +
  scale_x_break(c(-150,-1300))

top_100_vape_both %>% 
  ggplot(aes(x = log2FoldChange.12)) +
  geom_histogram(binwidth = 0.5, color = "white") +
  labs(x = "Estmate with Sample 12")

top_100_vape_both %>% 
  ggplot(aes(x = log2FoldChange.no12)) +
  geom_histogram(binwidth = 0.5, color = "white") +
  labs(x = "Estmate without Sample 12")

top_100_vape_both %>% 
  ggplot(aes(x = fold_change)) +
  geom_histogram(binwidth = 0.25, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "red") +
  labs(x = "Change in Estimate",
       y = "Count")

######################### Results Table log2 fold change ###############################
top_10_vape_p_change <- top_100_vape_both %>% 
  arrange(fold_percent_change) %>% 
  .[1:10,] %>% 
  select(row, log2FoldChange.12, log2FoldChange.no12, fold_percent_change)
top_10_vape_p_change <- left_join(top_10_vape_p_change, gene_annotations, by = c("row" = "gene"))
top_10_vape_p_change <- top_10_vape_p_change %>% 
  dplyr::mutate(log2FoldChange.12 = round(log2FoldChange.12, 3),
         log2FoldChange.no12 = round(log2FoldChange.no12, 3),
         fold_percent_change = format_num(fold_percent_change),
         Difference = format_num(log2FoldChange.12 - log2FoldChange.no12, digits = 3)) %>% 
  dplyr::select(row, symbol, log2FoldChange.12, log2FoldChange.no12, Difference, fold_percent_change) %>% 
  dplyr::rename("Ensemble ID" = row,
                "Gene Name" = symbol,
                "Log2 Fold Change (Model 1)" = log2FoldChange.12,
                "Log2 Fold Change (Model 2)" = log2FoldChange.no12,
                "% Change" = fold_percent_change)
  

#Write out results
# write_csv(top_10_vape_p_change, file = here("DataProcessed/sample12_sensitivity/top_10_genes_estchange_2022_05_10.csv"))
# write_csv(top_100_vape_both, file = here("DataProcessed/sample12_sensitivity/top_100_12sensitiv_2022_05_10.csv"))
# write_csv(vape_res, here("DataProcessed/sample12_sensitivity/Vape_Results_With_12_2022_05_11.csv"))
# write_csv(vape_res_no12, here("DataProcessed/sample12_sensitivity/Vape_Results_no12_2022_05_11.csv"))
