####################Load Libraries########################
library(tidyverse)
library(here)
library(diffEnrich)
library(clusterProfiler)

################## Collect and store pathways from KEGG ###################
kegg_hsa <- get_kegg('hsa', path = here("DataRaw/pathways"))

key <- kegg_hsa$ncbi_to_kegg
################## Read in tested genes################

################## Get Gene List and map to Entrez ID ####################
#Read it in
vape_res <- read_csv(here("DataProcessed/de_full/full_vape_res_2022_05_01.csv"))

#Apply cutoff
vape_res <- vape_res %>%
  filter(abs(log2FoldChange) > 2) %>% arrange(padj)

#Fix Ensemble ID's 
vape_res$ensg <- gsub("\\..*", "", vape_res$ensg)

#Get top 200
vape_top_200 <- vape_res[1:200,]

#Translate from Symbol to ENTREZID
ensembl_trans_200 <- bitr(vape_top_200$ensg, fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #Has only 1% fail to map
ensembl_trans <- bitr(vape_res$ensg, fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#Gene List
top_200_list<- ensembl_trans_200$ENTREZID
################# Subset gk obj #############################

ncbi_to_kegg <- kegg_hsa[["ncbi_to_kegg"]]
colnames(ncbi_to_kegg) <- c("ncbi_id", "gene")
ncbi_to_kegg$Entry <- gsub("ncbi-geneid:","",fixed = T, ncbi_to_kegg$ncbi_id)

ncbi_to_kegg_sub <- ncbi_to_kegg[ensembl_trans$ENTREZID %in% ncbi_to_kegg$Entry,] %>% 
  dplyr::select(-Entry)

colnames(ncbi_to_kegg_sub) <- c("V1", "V2")
kegg_hsa[["ncbi_to_kegg"]] <- ncbi_to_kegg_sub

enrich_try <- pathEnrich(gk_obj = kegg_hsa, gene_list = top_200_list, cutoff = 0.05, N = 2)

summary(enrich_try)

enrich_tab <- enrich_try$enrich_table