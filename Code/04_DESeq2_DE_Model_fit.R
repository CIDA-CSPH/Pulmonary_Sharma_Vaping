#Load libraries
library(ggpubr)
library(tidyverse)
library(here)
library(DESeq2)
library(pcaExplorer)
library(RColorBrewer)
library(kableExtra)

#read in files
filtered_gene_count <- as.data.frame(read_csv(here("DataProcessed/filtered_gene_count_2022_05_04.csv")))
genes <- scan(file = here("DataProcessed/filtered_gene_list_2022_04_20.txt"), character(), quote = "")
metadata_joined <- as.data.frame(read_csv(here("DataProcessed/metadata_joined_2022_04_20.csv")))
ruv_factor_dat <- read_csv(here("DataProcessed/ruv_factor_data_k2_2022_04_20.csv"))
ruv_counts<- read_csv(here("DataProcessed/RUV_k2_Counts_2022_05_01.csv"))
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

# #Try filter out Sample12
# #metadata
# metadata_joined <- metadata_joined %>%
#   filter(new_id != "Sample12")
# #count data
# filtered_gene_count <- filtered_gene_count %>%
#   select(-Sample12)

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

#Design Matrix
full_design <- ~vape_6mo_lab + sex_lab + age + recruitment_center + ruv_k1 + ruv_k2 

vape_only_reduced <- ~vape_6mo_lab + sex_lab + age + ruv_k1 + ruv_k2

center_only_reduced <- ~ sex_lab + age + recruitment_center + ruv_k1 + ruv_k2

reduced_design <- ~sex_lab + age + ruv_k1 + ruv_k2 

#Write Function to run DESeq
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

#Get the formatted results for p-value histograms

format_results <- function(deseq_run) {
  #Pull out resultsand sort by padj and pval
  deseq_res <- DESeq2::results(deseq_run, tidy = T, alpha = 0.05) %>% 
    arrange(padj, pvalue) %>% 
    tbl_df()
  return(deseq_res)
}

######################### Run DESeq ###############################
# Vape and Center
vape_center_de <- run_deseq_lrt(count_data = filtered_gene_count, 
                                 col_data = metadata_joined,
                                 full_mod = full_design,
                                 reduced_mod = reduced_design) 
#Vape Only
vape_de <- run_deseq_lrt(count_data = filtered_gene_count, 
                          col_data = metadata_joined,
                          full_mod = full_design,
                          reduced_mod = center_only_reduced)

#Center Only
center_de <- run_deseq_lrt(count_data = filtered_gene_count, 
                            col_data = metadata_joined,
                            full_mod = full_design,
                            reduced_mod = vape_only_reduced)

######################### Get Formatted Results ###############################
vape_center_res <- format_results(vape_center_de)

vape_res <- format_results(vape_de)

center_res <- format_results(center_de)

#Write out results
# write_csv(vape_center_res, file = here("DataProcessed/de_output/vape_center_res_2022_05_01.csv"))
# write_csv(vape_res, file = here("DataProcessed/de_output/vape_res_2022_05_01.csv"))
# write_csv(center_res, file = here("DataProcessed/de_output/center_res_2022_05_01.csv"))

######################### p-value histograms ###############################
#write a function
p_hist <- function(de_results) {
  temp_hist <- de_results %>% 
    ggplot(aes(x = pvalue))+
    geom_histogram()
  return(temp_hist)
}
#Vape and Center
vape_center_phist <- p_hist(vape_center_res) +
  labs(title = "Vape and Center")
#Vape Only
vape_only_phist <- p_hist(vape_res) +
  labs(title = "Vape Only")
#Center Only
center_only_phist <- p_hist(center_res) + 
  labs(title = "Center Only")

#fix metadata_joined vape status
metadata_joined$vape_6mo_lab <- if_else(metadata_joined$vape_6mo_lab == "vaped", "Vaped", "Not Vaped")
######################### Get top genes and tidy results for plotting with DE Norm Counts ###############################
gene_annotations <- gene_annotations %>% 
  mutate(gene = ENSG) %>% 
  select(gene, symbol)

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
#Vape and Center
vape_center_tcounts <- de_res_tidy(de_res = vape_center_res, 
                                   run_de = vape_center_de,
                                   col_data = metadata_joined)
#Vape only
vape_tcounts <- de_res_tidy(de_res = vape_res, 
                            run_de = vape_de,
                            col_data = metadata_joined)
#Center Only
center_tcounts <- de_res_tidy(de_res = center_res, 
                              run_de = center_de,
                              col_data = metadata_joined)

#write out files
# write_csv(vape_center_tcounts, file = here("DataProcessed/de_output/vape_center_de_normcounts_2022_05_01.csv"))
# write_csv(vape_tcounts, file = here("DataProcessed/de_output/vape_de_normcounts_2022_05_01.csv"))
# write_csv(center_tcounts, file = here("DataProcessed/de_output/center_de_normcounts_2022_05_01.csv"))

######################### Get top genes and tidy results for plotting with RUV Norm Counts ###############################
#fix ruv colnames
rownames(ruv_norm_counts) <- ruv_norm_counts$gene

ruv_norm_counts <- ruv_norm_counts %>% 
  select(-gene)

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

vape_center_ruv_tidy <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts,
                                de_res = vape_center_res,
                                col_data = metadata_joined)
vape_only_ruv_tidy <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts,
                                   de_res = vape_res,
                                   col_data = metadata_joined)

center_only_ruv_tidy <- ruv_res_tidy(ruv_count_dat = ruv_norm_counts,
                                     de_res = center_res,
                                     col_data = metadata_joined)

#write out results
# write_csv(vape_center_ruv_tidy, file = here("DataProcessed/de_output/vape_center_ruv_normcounts_2022_05_06.csv"))
# write_csv(vape_only_ruv_tidy, file = here("DataProcessed/de_output/vape_ruv_normcounts_2022_05_06.csv"))
# write_csv(center_only_ruv_tidy, file = here("DataProcessed/de_output/center_ruv_normcounts_2022_05_06.csv"))

######################### Top Gene Boxplots ###############################
#Fix labelling from metadata joined

#Set color pallettes to match
center_colors <- brewer.pal(10,"Paired")[4:6]
vape_colors <- brewer.pal(3,"Set2")
#Write a function
tcount_boxplot <- function(tcounts, mod) {
  #Stratified by Vape and Center
  vape_center_box <- tcounts %>% 
    ggplot(aes(vape_6mo_lab, expression, group = interaction(vape_6mo_lab, recruitment_center))) + 
    geom_boxplot(aes(fill = recruitment_center)) +
    geom_text(label = if_else(tcounts$Row.names == "Sample12", tcounts$Row.names, ""), col = "Red") +
    facet_wrap(~gene, scales="free_y") + 
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
    facet_wrap(~gene, scales="free_y") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom")
  
  #Stratified by center only
  center_box <- tcounts %>% 
    ggplot(aes(recruitment_center, expression, fill=recruitment_center)) + 
    geom_boxplot() +
    geom_text(label = if_else(tcounts$Row.names == "Sample12", tcounts$Row.names, ""), 
              col = "Red") +
    facet_wrap(~gene, scales="free_y") + 
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

#code to label sample 12
# geom_text(label = if_else(tcounts_vape_center$Row.names == "Sample12", tcounts_vape_center$Row.names, ""), col = "Red")

######################### Top Gene Boxplots DE Counts ###############################
#Vape and Center
vape_center_plots <- tcount_boxplot(tcounts = vape_center_tcounts, mod = "Vape and Center")

#Vape only
vape_plots <- tcount_boxplot(tcounts = vape_tcounts, mod = "Vape Only")

#Center Only
center_plots <- tcount_boxplot(tcounts = center_tcounts, mod = "Center Only")

#ruv
ruv_plots <- tcount_boxplot(tcounts = ruv_counts, mod = "RUV Only")

#Looking at plots
#Vape and Center
vape_center_plots[[1]]
vape_center_plots[[2]]
vape_center_plots[[3]]

#Vape only
vape_plots[[1]]
vape_plots[[2]]
vape_plots[[3]]

#Center Only
center_plots[[1]]
center_plots[[2]]
center_plots[[3]]

#RUV Only
ruv_plots[[1]]
ruv_plots[[2]]
ruv_plots[[3]]

######################### Top Gene Boxplots RUV Counts ###############################
#Vape and Center
vape_center_plots_ruv <- tcount_boxplot(tcounts = vape_center_ruv_tidy, mod = "Vape and Center")

#Vape only
vape_plots_ruv <- tcount_boxplot(tcounts = vape_only_ruv_tidy, mod = "Vape Only")

#Center Only
center_plots_ruv <- tcount_boxplot(tcounts = center_only_ruv_tidy, mod = "Center Only")


#Vape and Center
vape_center_plots_ruv[[1]]
vape_center_plots_ruv[[2]]
vape_center_plots_ruv[[3]]

#Vape only
vape_plots_ruv[[1]]
vape_plots_ruv[[2]]
vape_plots_ruv[[3]]

#Center Only
center_plots_ruv[[1]]
center_plots_ruv[[2]]
center_plots_ruv[[3]]

ggarrange(plotlist = list(vape_center_plots[[2]], vape_center_plots[[3]]), nrow = 1)



ggarrange(plotlist = ruv_plots, nrow = 3)

######################### Results Table ###############################
# function 
format_big <- function(x) {
  formatC(x, digits = 0, format = "f", big.mark = ",")
}

sig_gene_count <- function(de_res) {
  sum(de_res$padj < 0.05)
}

sig_genes_tab <- tibble(Model = c("Vape and Center", "Vape Only", "Center Only"),
                        "Significant Genes (p < 0.05)" = c(format_big(sig_gene_count(de_res = vape_center_res)),
                                                           format_big(sig_gene_count(de_res = vape_res)),
                                                           sig_gene_count(de_res = center_res))) %>% 
  kbl(digits = 3)



sum(sizeFactors(vape_center_de) == sizeFactors(vape_de))

counts(vape_center_de, normalized = T)[1,1]
counts(vape_de, normalized = T)[1,1]
