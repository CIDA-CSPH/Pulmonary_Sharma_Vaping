## ---------------------------
## Script name: DESEQ Results Formatting
##
## Purpose of script: Fomrat results for Sarah 
##
## Author: Trent Hawkins
##
## Date Created: 2022-12-01
## ---------------------------
## set working directory for Mac and PC

require(tidyverse)
require(here)
require(EnhancedVolcano)

results <- read_csv(here("DataProcessed/rna_seq/differential_expression/full_analysis/full_vape_res_2022_10_13.csv"))

results <- results %>%
  mutate(Signif = if_else(padj < 0.05, "*", "")) %>% 
  rename("ENSG" = ensg,
         "Symbol" = gene_name,
         "Log2FoldChange" = log2FoldChange,
         "p-value" = pvalue,
         "FDR" = padj) %>% 
  select(-c(gene_type, baseMean, lfcSE, stat)) %>% as.data.frame()

lab_italics <- paste0("italic('", results$Symbol, "')")
selectLab_italics = paste0(
  "italic('",
  c("NDC80", "PTGER3", "MUC5AC", "MUC12", "WNT5B", "WNT3A", "WNT4", 'IL1B', 'DPYS', 'DEFB4B', 'EIF1AY'),
  "')")

EnhancedVolcano(results,
                lab = lab_italics,
                x = 'Log2FoldChange',
                y = 'p-value',
                selectLab = selectLab_italics,
                xlab = bquote(~Log[2]~ 'FC'),
                ylab = bquote(~-Log[10]~ 'p-value'),
                pCutoff = 0.05,
                pCutoffCol = 'FDR',
                FCcutoff = 2.0,
                pointSize = 1.5,
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "",
                caption = "Total = 16,860 Genes")



write_csv(results, here("DataProcessed/rna_seq/differential_expression/full_analysis/vape_results_reporting_2022_12_01.csv"))
