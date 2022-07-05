####################Load Libraries########################
library(tidyverse)
library(here)
library(diffEnrich)
library(clusterProfiler)
library(fgsea)
library(reactome.db)
library(biomaRt)
library(kableExtra)

kablize <- function(tab, digits = 3) {
  kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("hold_position", "striped"), position = "center")
}

format_num <- function(number, digits = 0) {
  formatC(number, digits = digits, format = "f", big.mark = ",")
}
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
  dplyr::mutate(rank = sign(log2FoldChange)*-log10(pvalue)) %>%
  arrange(padj) %>% 
  dplyr::filter(!is.na(ENTREZID))

#ENTRES
ranks_entrez <- vape_res$rank
names(ranks_entrez) <- vape_res$ENTREZID

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

kegg_key <- left_join(vape_res, ncbi_to_name, by = c("ENTREZID" = "ncbi"))

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
# hsa_ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

go_key <- read_rds(here("DataRaw/pathways/go_mart_2022_06_20.rds"))


# go_key <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name",
#                                "go_id", "namespace_1003", "name_1006", "definition_1006"),
#                 filters = "ensembl_gene_id",
#                 values = vape_res$ensg,
#                 mart = hsa_ensembl)
# write_rds(go_key, file = here("DataRaw/pathways/go_mart_2022_06_20.rds"))

go_key <- go_key %>% 
  filter(!is.na(entrezgene_id),
         name_1006 != "")

go_key_mapped <- go_key %>% 
  dplyr::group_by(name_1006) %>% 
  dplyr::summarise(NCBI = toString(entrezgene_id),
                   Gene_Name = toString(external_gene_name)) %>% 
  dplyr::ungroup()

go_pathways <- as.list(strsplit(go_key_mapped$NCBI, split = ", "))

names(go_pathways) <- go_key_mapped$name_1006

go_pathways <- go_pathways[!duplicated(names(go_pathways))]

########################## GSEA ####################################
#Run GSEA
set.seed(404)

#Kegg
kegg_res <- fgsea(pathways = kegg_pathways,
                    stats = ranks_entrez,
                    minSize = 2,
                    maxSize = 500,
                    nPermSimple = 100000)
#React
react_res <- fgsea(pathways = react_path,
                     stats = ranks_entrez,
                     minSize = 2,
                     maxSize = 500,
                     nPermSimple = 100000)
#GO
go_res <- fgsea(pathways = go_pathways,
                  stats = ranks_entrez,
                  minSize = 2,
                  maxSize = 500,
                  nPermSimple = 100000)

res_collapse <- collapsePathways(fgseaRes = kegg_res, 
                                 pathways = kegg_pathways,
                                 stats = ranks_entrez)
###################### Match up/clean/format results with Gene Names############

format_res <- function(result, path, ranks) {
  set.seed(404)
  res_collapse <- collapsePathways(fgseaRes = result, 
                                   pathways = path,
                                   stats = ranks)
  
  res_ind <- result[pathway %in% res_collapse$mainPathways,]
  
  entrez_genes <- unique(unlist(res_ind[["leadingEdge"]]))
  
  symbol_map <- bitr(entrez_genes, fromType = "ENTREZID", 
                     toType = "SYMBOL", OrgDb="org.Hs.eg.db")
  
  positive_entrez <- NULL
  
  negative_entrez <- NULL
  
  positive_symbol <- NULL
  
  negative_symbol <- NULL
  
  positive_length <- NULL
  
  negative_length <- NULL
  
  for (i in 1:nrow(res_ind)) {
    entrez_temp <- unique(unlist(res_ind[i,"leadingEdge"]))
    
    vape_res_temp <- vape_res[vape_res$ENTREZID %in% entrez_temp,]
    
    positive_genes_temp <- vape_res_temp[vape_res_temp$log2FoldChange > 0, c('ENTREZID','gene_name')]
    
    negative_genes_temp <- vape_res_temp[vape_res_temp$log2FoldChange < 0, c('ENTREZID','gene_name')]
    
    positive_entrez <- append(positive_entrez, list(positive_genes_temp$ENTREZID))
    
    negative_entrez <- append(negative_entrez, list(negative_genes_temp$ENTREZID))
    
    positive_symbol <- append(positive_symbol, list(positive_genes_temp$gene_name))
    
    negative_symbol <- append(negative_symbol, list(negative_genes_temp$gene_name))
    
    positive_length <- append(positive_length, length(positive_genes_temp$ENTREZID))
    
    negative_length <- append(negative_length, length(negative_genes_temp$ENTREZID))
  
    print(positive_length)
  }
  
  res_ind[,'positive_entrez'] <- positive_entrez

  res_ind[,'positive_symbol'] <- positive_symbol

  res_ind[,'negative_entrez'] <- negative_entrez

  res_ind[,'negative_symbol'] <- negative_symbol
  
  res_ind[,'positive_length'] <- positive_length
  
  res_ind[,'negative_length'] <- negative_length

  res_ind <- res_ind %>%
    dplyr::select(pathway, pval, padj, ES, NES, size, leadingEdge, positive_entrez, positive_symbol,positive_length,
                  negative_entrez, negative_symbol, negative_length)
  return(res_ind)}


kegg_ind <- format_res(result = kegg_res,
                       path = kegg_pathways,
                       ranks = ranks_entrez)

react_ind <- format_res(result = react_res,
                        path = react_path,
                        ranks = ranks_entrez)

go_ind <- format_res(result = go_res,
                     path = go_pathways,
                     ranks = ranks_entrez)


########################## Visualizations and Tables ###########################

get_paths_res <- function(res, res_ind) {
  paths_tested <- as.numeric(nrow(res))
  
  ind_paths <- nrow(res_ind)
  
  genes_matched <- length(unique(unlist(res[["leadingEdge"]])))
  
  sig_paths_05 <- as.numeric(nrow(res_ind %>% dplyr::filter(padj < 0.05)))
  
  sig_paths_10 <- as.numeric(nrow(res_ind %>% dplyr::filter(padj < 0.1)))
  
  sig_paths_20 <- as.numeric(nrow(res_ind %>% dplyr::filter(padj < 0.2)))
  
  return(c(paths_tested, genes_matched, ind_paths, sig_paths_05, sig_paths_10, sig_paths_20))
}

kegg_sum <- get_paths_res(kegg_res, kegg_ind)

react_sum <- get_paths_res(react_res, react_ind)

go_sum <- get_paths_res(go_res, go_ind)

paths_res <- rbind(kegg_sum, react_sum, go_sum)

rownames(paths_res) <- c("KEGG", "Reactome", "GO")
colnames(paths_res) <- c("Paths Tested",  "Genes Matched", "Independent Pathways", "FDR < 0.05", "FDR < 0.1", "FDR < 0.2")

#Summary Stats Table
paths_res %>% 
  kablize()

kegg_ind$pathway <- gsub(".{23}$", "", kegg_ind$pathway)

kegg_res$pathway <- gsub(".{23}$", "", kegg_res$pathway)


#Top Results for each pathway

format_tab <- function(res_ind) {
  res_ind %>% 
    arrange(padj, pval) %>%
    dplyr::select(-c(leadingEdge, positive_entrez, positive_symbol, negative_entrez, negative_symbol)) %>% 
    dplyr::mutate(ES = round(ES, 2),
                  NES = round(NES, 2),
                  pval = if_else(pval < 0.001, '<0.001', as.character(round(pval, 3))),
                  padj = if_else(padj < 0.001, '<0.001', as.character(round(padj, 3)))) %>%
    dplyr::rename(Pathway = pathway,
                  FDR = padj,
                  "Up-Regulated Genes" = positive_length,
                  "Down-Regulated Genes" = negative_length) %>% 
    
    .[1:10] %>% 
    kablize()
}

format_tab(kegg_ind)

format_tab(react_ind)

format_tab(go_ind)




plotEnrichment(kegg_pathways, ranks_entrez) +
  labs(title = "KEGG Enrichment Score Across Ranks")

plotEnrichment(react_path, ranks_entrez) +
  labs(title = "Reactome Enrichment Score Across Ranks")

plotEnrichment(go_pathways, ranks_entrez) +
  labs(title = "GO Enrichment Score Across Ranks")

plotGseaTable(kegg_pathways[kegg_ind$pathway[1:10]], ranks_entrez, kegg_ind)

plotGseaTable(react_path[react_ind$pathway[1:10]], ranks_entrez, react_ind)

plotGseaTable(go_pathways[go_ind$pathway[1:10]], ranks_entrez, go_ind)
