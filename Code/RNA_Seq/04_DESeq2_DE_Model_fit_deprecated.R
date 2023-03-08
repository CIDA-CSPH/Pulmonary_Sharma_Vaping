######################### Load Libraries ###############################
library(tidyverse)
library(here)
library(DESeq2)
library(RColorBrewer)
library(kableExtra)

######################### Read in files so that they are unaltered and join ###############################
filtered_gene_count <- read.csv(here("DataProcessed/rna_seq/differential_expression/full_analysis/filtered_gene_count_2022_10_13.csv"), 
                                row.names = 1, header = T) 
metadata_joined <- read.csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv"))
ruv_factor_dat <- read.csv(here("DataProcessed/rna_seq/ruv/ruv_factor_data_k2_2022_10_13.csv")) %>% 
  select(c(rna_id, W_1, W_2))
ruv_norm_counts <- read.csv(here("DataProcessed/rna_seq/ruv/RUV_k2_norm_counts_2022_10_13.csv"), 
                            row.names = 1, header = T)
gene_annotations <- read_tsv(here("DataRaw/RNA_Seq/gencode_annotations_choo.txt"))

#only metadata with rna_seq dat
metadata_joined <- metadata_joined %>% 
  filter(rna_id %in% names(filtered_gene_count))
  
#join metadata with ruv factor data
metadata_joined <- left_join(metadata_joined, ruv_factor_dat, by = "rna_id") %>% 
  dplyr::rename("ruv_k1" = W_1,
                "ruv_k2" = W_2) 

rownames(metadata_joined) <- metadata_joined$rna_id

######################### Prepare Metadata ###############################
#Set factor levels
metadata_joined$vape_6mo_lab <- factor(metadata_joined$vape_6mo_lab, levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"), 
                                       labels = c("not_vaped", "vaped"))

metadata_joined$sex_lab <- factor(metadata_joined$sex_lab, levels = c("Female", "Male"))

metadata_joined$recruitment_center <- factor(metadata_joined$recruitment_center, levels = c("Aurora", "CommCity/Denver", "Pueblo"), 
                                             labels = c("Aurora", "CommCity_Denver", "Pueblo"))
#Scaling Age
metadata_joined$age <- scale(metadata_joined$age)

######################### Design Matrices ###############################
full_design <- ~ sex_lab + age + recruitment_center + ruv_k1 + ruv_k2 + vape_6mo_lab 

test_center <- ~ sex_lab + age + ruv_k1 + ruv_k2 + vape_6mo_lab

test_vape <- ~ sex_lab + age + recruitment_center + ruv_k1 + ruv_k2

test_vape_center <- ~sex_lab + age + ruv_k1 + ruv_k2 
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
  deseq_res <- DESeq2::results(deseq_run, tidy = T, alpha = 0.05, name = "vape_6mo_lab_vaped_vs_not_vaped") %>% 
    arrange(padj, pvalue) %>% 
    as_tibble() %>% 
    left_join(.,annotation, by = c('row' = 'ENSG')) %>% 
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
  genes_of_interest <- de_res$ensg[1:100]
  
  #Normalize Counts
  run_de <- estimateSizeFactors(run_de)
  
  #Join and tidy
  norm_counts<- t(log10((counts(run_de[genes_of_interest, ], 
                                normalized=TRUE)+.5))) %>%
    merge(col_data, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(genes_of_interest) + 1):ncol(.))
  
  #add symbol name
  gene_symbols <- gene_annotations[gene_annotations$ENSG %in% genes_of_interest,]
  
  norm_counts <- left_join(norm_counts, gene_symbols, by = c("gene" = "ENSG")) 
  
  return(norm_counts)
}

#tidy results for plotting (RUV-Normalized)
ruv_res_tidy <- function(ruv_count_dat, de_res, col_data, annotation) {
  #Get top 4 genes of interest
  genes_of_interest <- de_res$ensg[1:100]
  
  #Join and tidy
  norm_counts<- t(log10((ruv_count_dat[genes_of_interest,] +.5))) %>%
    merge(col_data, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-length(genes_of_interest) + 1):ncol(.))
  
  print(colnames(norm_counts))
  print(colnames(annotation))
  #join back with gene name
  norm_counts <- left_join(norm_counts, annotation, by = c("gene" = "ENSG"))
  
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
    scale_fill_manual(values = center_colors)
  
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
vape_center_de <- run_deseq_lrt(count_data = filtered_gene_count, 
                                col_data = metadata_joined,
                                full_mod = full_design,
                                reduced_mod = test_vape_center) 
#Vape Only
vape_de <- run_deseq_lrt(count_data = filtered_gene_count, 
                         col_data = metadata_joined,
                         full_mod = full_design,
                         reduced_mod = test_vape)

#Center Only
center_de <- run_deseq_lrt(count_data = filtered_gene_count, 
                           col_data = metadata_joined,
                           full_mod = full_design,
                           reduced_mod = test_center)

######################### Get Formatted Results ###############################
vape_center_res <- format_results(vape_center_de, gene_annotations)

vape_res <- format_results(vape_de, gene_annotations)

center_res <- format_results(center_de, gene_annotations)

#Write out results
write_csv(vape_center_res, file = here("DataProcessed/rna_seq/differential_expression/full_analysis/vape_center_res_2022_10_13.csv"))
write_csv(vape_res, file = here("DataProcessed/rna_seq/differential_expression/full_analysis/full_vape_res_2022_10_13.csv"))
write_csv(center_res, file = here("DataProcessed/rna_seq/differential_expression/full_analysis/center_res_2022_10_13.csv"))

######################### P-Value Histograms ###############################
#Vape and Center
vape_center_phist <- p_hist(vape_center_res) +
  labs(title = "Vape and Center")
#Vape Only
vape_only_phist <- p_hist(vape_res) +
  labs(title = "Vape Only")
#Center Only
center_only_phist <- p_hist(center_res) + 
  labs(title = "Center Only")

######################### Top Genes Tidy DE-Norm Counts ###############################
#fix metadata_joined vape status
metadata_joined$vape_6mo_lab <- if_else(metadata_joined$vape_6mo_lab == "vaped", "Vaped", "Not Vaped")

#Vape and Center
vape_center_de_tidy <- de_res_tidy(de_res = vape_center_res, 
                                   run_de = vape_center_de,
                                   col_data = metadata_joined)
#Vape only
vape_de_tidy  <- de_res_tidy(de_res = vape_res, 
                             run_de = vape_de,
                             col_data = metadata_joined)
#Center Only
center_de_tidy  <- de_res_tidy(de_res = center_res, 
                               run_de = center_de,
                               col_data = metadata_joined)

#write out files
# write_csv(vape_center_tcounts, file = here("DataProcessed/de_output/vape_center_de_normcounts_2022_05_01.csv"))
# write_csv(vape_tcounts, file = here("DataProcessed/de_output/vape_de_normcounts_2022_05_01.csv"))
# write_csv(center_tcounts, file = here("DataProcessed/de_output/center_de_normcounts_2022_05_01.csv"))

######################### Top Genes Tidy RUV-Norm Counts ###############################
vape_center_ruv_tidy <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts,
                                     de_res = vape_center_res,
                                     col_data = metadata_joined, 
                                     annotation = gene_annotations)

vape_ruv_tidy <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts,
                                   de_res = vape_res,
                                   col_data = metadata_joined,
                                   annotation = gene_annotations)


center_ruv_tidy <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts,
                                de_res = center_res,
                                col_data = metadata_joined,
                                annotation = gene_annotations)

#write out results
write_csv(vape_center_ruv_tidy, file = here("DataProcessed/rna_seq/differential_expression/full_analysis/vape_center_ruv_normcounts_2022_10_13.csv"))
write_csv(vape_ruv_tidy, file = here("DataProcessed/rna_seq/differential_expression/full_analysis/vape_ruv_normcounts_2022_10_13.csv"))
write_csv(center_ruv_tidy, file = here("DataProcessed/rna_seq/differential_expression/full_analysis/center_ruv_normcounts_2022_10_13.csv"))

######################### Top Gene Boxplots DE Counts ###############################
#Set color pallettes to match
center_colors <- brewer.pal(10,"Paired")[4:6]
vape_colors <- brewer.pal(3,"Set2")

#Vape and Center
vape_center_plots <- tcount_boxplot(tcounts = vape_center_de_tidy, mod = "Vape and Center")

#Vape only
vape_plots <- tcount_boxplot(tcounts = vape_de_tidy, mod = "Vape Only")

#Center Only
center_plots <- tcount_boxplot(tcounts = center_de_tidy, mod = "Center Only")

#Looking at plots
#Vape and Center
vape_center_plots[[1]]
vape_plots[[1]]
center_plots[[1]]

#Vape only
vape_center_plots[[2]]
vape_plots[[2]]
center_plots[[2]]

#Center Only
vape_center_plots[[3]]
vape_plots[[3]]
center_plots[[3]]

######################### Top Gene Boxplots RUV Counts ###############################
#Vape and Center
vape_center_plots_ruv <- tcount_boxplot(tcounts = vape_center_ruv_tidy, mod = "Vape and Center")

#Vape only
vape_plots_ruv <- tcount_boxplot(tcounts = vape_ruv_tidy, mod = "Vape Only")

#Center Only
center_plots_ruv <- tcount_boxplot(tcounts = center_ruv_tidy, mod = "Center Only")


#Vape and Center
vape_center_plots_ruv[[1]]
vape_plots_ruv[[1]]
center_plots_ruv[[1]]

#Vape only
vape_center_plots_ruv[[2]]
vape_plots_ruv[[2]]
center_plots_ruv[[2]]

#Center Only
vape_center_plots_ruv[[3]]
vape_plots_ruv[[3]]
center_plots_ruv[[3]]

######################### Results Table ###############################
sig_genes_tab <- tibble(Model = c("Vape and Center", "Vape Only", "Center Only"),
                        "Significant Genes (p < 0.05)" = c(format_num(sig_gene_count(de_res = vape_center_res)),
                                                           format_num(sig_gene_count(de_res = vape_res)),
                                                           sig_gene_count(de_res = center_res))) %>% 
  kbl(digits = 3)


