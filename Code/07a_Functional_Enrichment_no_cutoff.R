####################Load Libraries########################
library(tidyverse)
library(here)
library(diffEnrich)
library(clusterProfiler)

################## Collect and store pathways from KEGG ###################
#kegg_hsa <- get_kegg('hsa', path = here("DataRaw/pathways"))

ncbi_to_kegg <- read_tsv("DataRaw/pathways/ncbi_to_kegg2022-06-15Release_102.0+_06-15_Jun_22.txt")
kegg_to_pathway <- read_tsv("DataRaw/pathways/kegg_to_pathway2022-06-15Release_102.0+_06-15_Jun_22.txt")
pathway_to_species <- read_tsv("DataRaw/pathways/pathway_to_species2022-06-15Release_102.0+_06-15_Jun_22.txt")

kegg_hsa <- list(ncbi_to_kegg, kegg_to_pathway, pathway_to_species)
names(kegg_hsa) <- c("ncbi_to_kegg", "kegg_to_pathway", "pathway_to_species")
################## Read in tested genes################

################## Get Gene List and map to Entrez ID ####################
##### FULL LIST
#Read it in
vape_res <- read_csv(here("DataProcessed/de_full/full_vape_res_2022_05_01.csv"))

#Fix Ensemble ID's 
vape_res$ensg <- gsub("\\..*", "", vape_res$ensg)

#translate to ENTREZ
ensembl_trans <- bitr(vape_res$ensg, fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

#### TAKE TOP 200
#Get top 200
vape_top_200 <- vape_res[1:200,]

#Translate from Symbol to ENTREZID
ensembl_trans_200 <- bitr(vape_top_200$ensg, fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #Has only 1% fail to map

#Gene List
top_200_list<- ensembl_trans_200$ENTREZID
################# Subset gk obj #############################

ncbi_to_kegg <- kegg_hsa[["ncbi_to_kegg"]]
colnames(ncbi_to_kegg) <- c("ncbi_id", "gene")
ncbi_to_kegg$Entry <- gsub("ncbi-geneid:","",fixed = T, ncbi_to_kegg$ncbi_id)

ncbi_to_kegg_sub <- ncbi_to_kegg[ncbi_to_kegg$Entry %in% ensembl_trans$ENTREZID,] %>% 
  dplyr::select(-Entry)

colnames(ncbi_to_kegg_sub) <- c("V1", "V2")
kegg_hsa[["ncbi_to_kegg"]] <- ncbi_to_kegg_sub

enrich_try <- pathEnrich(gk_obj = kegg_hsa, gene_list = top_200_list, cutoff = 0.05, N = 2)

summary(enrich_try)

enrich_tab <- enrich_try$enrich_table