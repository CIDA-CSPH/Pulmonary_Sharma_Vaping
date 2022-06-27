####################Load Libraries########################
library(tidyverse)
library(here)
library(diffEnrich)
library(clusterProfiler)
library(fgsea)
library(reactome.db)
library(biomaRt)

################## Read in results #################
vape_res <- read_csv(here("DataProcessed/de_full/full_vape_res_2022_05_01.csv"))

#Fix Ensemble ID's  and join with symbols
vape_res$ensg <- gsub("\\..*", "", vape_res$ensg)

#translate to ENTREZ
genes_tested <- bitr(vape_res$ensg, fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

vape_res <- left_join(vape_res, genes_tested, by = c("ensg" = "ENSEMBL"))

################# log2(FC) > 2 ############################
vape_res <- vape_res %>% 
  dplyr::filter(abs(log2FoldChange) > 2) %>% 
  dplyr::arrange(padj)

################# Get Ranks ############################
#Make Ranks
vape_res <- vape_res %>% 
  dplyr::mutate(rank = sign(log2FoldChange)*-log10(padj)) %>%
                  arrange(padj)

vape_res_entrez_only <- vape_res %>% 
  dplyr::filter(!is.na(ENTREZID))

#Subset to named vector
#ENSEMBL
ranks <- vape_res$rank
names(ranks) <- vape_res$ensg

#ENTRES
ranks_entrez <- vape_res_entrez_only$rank
names(ranks_entrez) <- vape_res_entrez_only$ENTREZID

############################### KEGG ###########################################
#get_kegg("hsa", path = here("DataRaw/pathways"))

ncbi_to_kegg <- read_tsv(here("DataRaw/pathways/ncbi_to_kegg2022-06-20Release_102.0+_06-20_Jun_22.txt"), col_names = F)
kegg_to_path <- read_tsv(here("DataRaw/pathways/kegg_to_pathway2022-06-20Release_102.0+_06-20_Jun_22.txt"), col_names = F)
path_to_name <- read_tsv(here("DataRaw/pathways/pathway_to_species2022-06-20Release_102.0+_06-20_Jun_22.txt"), col_names = F)

#set column names
colnames(ncbi_to_kegg) <- c("ncbi", "kegg")
colnames(kegg_to_path) <- c("kegg", "path")
colnames(path_to_name) <- c("path", "name")


ncbi_to_path <- left_join(ncbi_to_kegg, kegg_to_path, by = "kegg")
ncbi_to_name <- left_join(ncbi_to_path, path_to_name, by = "path")

#Produce Kegg Key
ncbi_to_name <- ncbi_to_name %>% 
  dplyr::mutate(ncbi = gsub("ncbi-geneid:", "", ncbi),
                kegg = gsub("hsa:", "", kegg),
                path = gsub("path:", "", path))

kegg_key <- left_join(vape_res_entrez_only, ncbi_to_name, by = c("ENTREZID" = "ncbi"))

kegg_key <- kegg_key %>% 
  dplyr::select(ensg, ENTREZID, gene_name, path, name)

#Extract pathways and map them 
pathways <- as.list(unique(kegg_key$name[!is.na(kegg_key$name)]))

kegg_key_mapped <- kegg_key %>% 
  dplyr::filter(!is.na(name)) %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(ENSG = toString(ensg),
                   NCBI = toString(ENTREZID),
                   Gene_Name = toString(gene_name)) %>% 
  dplyr::ungroup()

kegg_pathways <- as.list(strsplit(kegg_key_mapped$NCBI, split = ", "))

names(kegg_pathways) <- kegg_key_mapped$name
############################ Reactome ###############################
react_path <- reactomePathways(names(ranks_entrez))

############################# Gene Ontology ####################################
#hsa_ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

go_key <- read_rds(here("DataRaw/pathways/go_mart_2022_06_20.rds"))


# go_key <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name",
#                                "go_id", "namespace_1003", "name_1006", "definition_1006"),
#                 filters = "ensembl_gene_id",
#                 values = vape_res$ensg,
#                 mart = hsa_ensembl)
# write_rds(go_key, file = here("DataRaw/pathways/go_mart_2022_06_20.rds"))

go_key <- go_key %>% 
  filter(name_1006 != "")

go_key_mapped <- go_key %>% 
  dplyr::group_by(name_1006) %>% 
  dplyr::summarise(ENSG = toString(ensembl_gene_id),
                   NCBI = toString(entrezgene_id),
                   Gene_Name = toString(external_gene_name)) %>% 
  dplyr::ungroup()

go_pathways <- as.list(strsplit(go_key_mapped$NCBI, split = ", "))

names(go_pathways) <- go_key_mapped$name_1006

go_pathways <- go_pathways[!duplicated(names(go_pathways))]
########################## GSEA ####################################
run_gsea <- function(top_genes, permutations, ranks) {
  set.seed(404)
  
  kegg_res <- fgsea(pathways = kegg_pathways,
                    stats = ranks,
                    minSize = 2,
                    maxSize = 500,
                    nPermSimple = 10000)
  react_res <- fgsea(pathways = react_path,
                     stats = ranks,
                     minSize = 2,
                     maxSize = 500,
                     nPermSimple = 100000)
  go_res <- fgsea(pathways = go_pathways,
                  stats = ranks,
                  minSize = 2,
                  maxSize = 500,
                  nPermSimple = 100000)
}


x <- list(kegg_pathways, react_path)
x <- base::append(x, list(go_pathways))

#write_rds(kegg_res, here("DataProcessed/gsea_res/kegg_res_ENTREZ_2022_06_21.rds"))
#write_rds(react_res, here("DataProcessed/gsea_res/react_res_ENTREZ_2022_06_21.rds"))
#write_rds(go_res, here("DataProcessed/gsea_res/go_res_ENTREZ_2022_06_21.rds"))

kegg_res_ensg <- read_rds(here("DataProcessed/gsea_res/kegg_res_2022_06_21.rds"))
go_res_ensg <- read_rds(here("DataProcessed/gsea_res/go_res_2022_06_21.rds"))

kegg_res_entrez <- read_rds(here("DataProcessed/gsea_res/kegg_res_ENTREZ_2022_06_21.rds"))
react_res_entrez <- read_rds(here("DataProcessed/gsea_res/react_res_ENTREZ_2022_06_21.rds"))
go_res_entrez <- read_rds(here("DataProcessed/gsea_res/go_res_ENTREZ_2022_06_21.rds"))

#Combine all entrez results
kegg_res_entrez$db <- "KEGG"
react_res_entrez$db <- "REACT"
go_res_entrez$db <- "GO"

all_res_entrez <- rbind(kegg_res_entrez, react_res_entrez, go_res_entrez)
########################## Visualizations and Tables ###########################
#Find Genes tested for each model
genes_tested_Entrez <- as.numeric(nrow(vape_res_entrez_only))

genes_tested_ensg <- as.numeric(nrow(vape_res))

#Paths tested and significant paths for each model
paths_res_entrez <- matrix(nrow = 3, ncol = 4) 
rownames(paths_res_entrez) <- c("KEGG", "Reactome", "GO")
colnames(paths_res_entrez) <- c("Paths_Tested", "Sig_Paths_05", "Sig_Paths_10", "Sig_Paths_20")

rownum <- 1

for (i in list(kegg_res_entrez, react_res_entrez, go_res_entrez)) {
  paths_tested <- as.numeric(nrow(i))
  
  sig_paths_05 <- as.numeric(nrow(i %>% dplyr::filter(padj < 0.05)))
  
  sig_paths_10 <- as.numeric(nrow(i %>% dplyr::filter(padj < 0.1)))
  
  sig_paths_20 <- as.numeric(nrow(i %>% dplyr::filter(padj < 0.2)))
  
  paths_res_entrez[rownum,] <- c(paths_tested, sig_paths_05, sig_paths_10, sig_paths_20)
  
  rownum <- rownum + 1
}


                                     
                                     