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
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}

format_num <- function(number, digits = 0, format = "f") {
  formatC(number, digits = digits, format = format, big.mark = ",")
}
################## Read in results #################
vape_res <- read_csv(here("DataProcessed/rna_seq/differential_expression/full_analysis/archive/full_vape_res_2022_10_06.csv"))

#Fix Ensemble ID's  and join with symbols
vape_res$ensg <- gsub("\\..*", "", vape_res$ensg)

################# log2(FC) > 2 ############################
vape_res <- vape_res %>% 
  dplyr::filter(abs(log2FoldChange) > 2) %>% 
  dplyr::arrange(padj)

#translate to ENTREZ
genes_tested <- bitr(vape_res$ensg, fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

vape_res <- left_join(vape_res, genes_tested, by = c("ensg" = "ENSEMBL"))



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
#get_kegg("hsa", path = here("DataRaw/enrichment"))

ncbi_to_kegg <- read_tsv("DataRaw/enrichment/ncbi_to_kegg2022-06-20Release_102.0+_06-20_Jun_22.txt")
kegg_to_path <- read_tsv("DataRaw/enrichment/kegg_to_pathway2022-06-20Release_102.0+_06-20_Jun_22.txt")
path_to_name <- read_tsv("DataRaw/enrichment/pathway_to_species2022-06-20Release_102.0+_06-20_Jun_22.txt")

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
#hsa_ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

go_key <- read_rds(here("DataRaw/enrichment/go_mart_2022_06_20.rds"))


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
                    nPermSimple = 1000)
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



###################### Match up/clean/format results with Gene Names############

format_res <- function(result, path, ranks) {
  set.seed(404)
  res_collapse <- collapsePathways(fgseaRes = result, 
                                   pathways = path,
                                   stats = ranks)
  
  #Save the independent path that the dependent path collapses into
  ind_path <- tibble(pathway = names(res_collapse$parentPathways),
                     independent_path = res_collapse$parentPathways)
  
  #map to the independent pathway
  result <- left_join(result, ind_path, by = "pathway")
  
  #Map entrez ID's to symbol
  entrez_genes <- unique(unlist(result[["leadingEdge"]]))
  
  symbol_map <- bitr(entrez_genes, fromType = "ENTREZID", 
                     toType = "SYMBOL", OrgDb="org.Hs.eg.db")
  
  #Prepare vectors to store gene info
  positive_entrez <- NULL
  
  negative_entrez <- NULL
  
  positive_symbol <- NULL
  
  negative_symbol <- NULL
  
  positive_length <- NULL
  
  negative_length <- NULL
  
  leadingEdge_unlist <- NULL
  
  leadingEdge_symbol <- NULL
  
  for (i in 1:nrow(result)) {
    #Get the entrez gene list for the pathway
    entrez_temp <- unique(unlist(result[i,"leadingEdge"]))
    #Subset the diff exp results 
    vape_res_temp <- vape_res[vape_res$ENTREZID %in% entrez_temp,]
    
    #All genes_temp
    all_genes_temp <- vape_res_temp[,c('ENTREZID','gene_name')]
    
    #Figure out if log2 FC is + or -
    positive_genes_temp <- vape_res_temp[vape_res_temp$log2FoldChange > 0, c('ENTREZID','gene_name')]
    
    negative_genes_temp <- vape_res_temp[vape_res_temp$log2FoldChange < 0, c('ENTREZID','gene_name')]
    
    #All Genes
    leadingEdge_symbol <- append(leadingEdge_symbol, paste0(all_genes_temp$gene_name, collapse = ', '))
    
    #Positive Genes
    positive_entrez <- append(positive_entrez, paste0(positive_genes_temp$ENTREZID, collapse = ', '))
    
    positive_symbol <- append(positive_symbol, paste0(positive_genes_temp$gene_name, collapse = ', '))
    
    positive_length <- append(positive_length, length(positive_genes_temp$ENTREZID))
    
    #Negative Genes
    negative_entrez <- append(negative_entrez, paste0(negative_genes_temp$ENTREZID, collapse = ', '))
  
    negative_symbol <- append(negative_symbol, paste0(negative_genes_temp$gene_name, collapse = ', '))
    
    negative_length <- append(negative_length, length(negative_genes_temp$ENTREZID))
    
    #fix leading edge
    leadingEdge_unlist <- append(leadingEdge_unlist, paste0(unlist(result$leadingEdge[i]), collapse = ', '))
  }
  
  #Add categorized genes to result
  result[,'positive_entrez'] <- positive_entrez

  result[,'positive_symbol'] <- positive_symbol

  result[,'negative_entrez'] <- negative_entrez

  result[,'negative_symbol'] <- negative_symbol
  
  result[,'positive_length'] <- positive_length
  
  result[,'negative_length'] <- negative_length
  
  result[,'leadingEdge'] <- leadingEdge_unlist
  
  result[,'leadingEdge_symbol'] <- leadingEdge_symbol

  result <- result %>%
    dplyr::select(pathway, independent_path, pval, padj, ES, NES, size, leadingEdge, leadingEdge_symbol,
                  positive_entrez, positive_symbol, positive_length,
                  negative_entrez, negative_symbol, negative_length)
  return(result)}




kegg_ind <- format_res(result = kegg_res,
                       path = kegg_pathways,
                       ranks = ranks_entrez)

react_ind <- format_res(result = react_res,
                        path = react_path,
                        ranks = ranks_entrez)

go_ind <- format_res(result = go_res,
                     path = go_pathways,
                     ranks = ranks_entrez)

write_csv(kegg_ind, here("DataProcessed/gsea/kegg_results_2022_10_06.csv"))
write_csv(react_ind, here("DataProcessed/gsea/react_results_2022_10_06.csv"))
write_csv(go_ind, here("DataProcessed/gsea/go_results_2022_10_06.csv"))
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

kegg_ind$independent_path <- gsub(".{23}$", "", kegg_ind$independent_path)



#Top Results for each pathway

format_tab <- function(res_ind) {
  
  res_ind <- res_ind %>% 
    arrange(padj, pval)
  
  res_ind %>% 
    dplyr::select(-c(leadingEdge, positive_entrez, positive_symbol, negative_entrez, negative_symbol)) %>% 
    dplyr::mutate(ES = round(ES, 2),
                  NES = round(NES, 2),
                  pval = if_else(pval < 0.001, format_num(pval, digits = 2, format = "e"), as.character(round(pval, 5))),
                  padj = if_else(padj < 0.001, format_num(padj, digits = 2, format = "e"), as.character(round(padj, 5))),
                  independent_path = if_else(is.na(independent_path), "Ind", independent_path)) %>%
    dplyr::rename(Pathway = pathway,
                  "Independent Pathway" = independent_path,
                  FDR = padj,
                  "Up-Regulated Genes" = positive_length,
                  "Down-Regulated Genes" = negative_length) 
}

kegg_res_supplement <- format_tab(kegg_ind) 

react_res_supplement <- format_tab(react_ind)

go_res_suplement <- format_tab(go_ind)

################# Format results for supplementary table ######################
sup_table_maker <- function(res) {
  #Make a row for each leading edge gene
  long_tab <- res %>% 
    separate_rows(leadingEdge, sep = ', ') %>% 
    dplyr::arrange(padj)
  
  #For each gene find if it is up- or down-regulated
  regulation <- NULL
  
  vape_res_rows <- match(as.numeric(long_tab$leadingEdge),
                               vape_res$ENTREZID)
  
  symbol_match <- vape_res[vape_res_rows, "gene_name"]
  
  
  #Format
  long_tab <- long_tab %>% 
    dplyr::select(-c(9:14)) %>% 
    dplyr::mutate(regulation = if_else(vape_res[vape_res_rows, "log2FoldChange"] > 0, '+', '-'),
                  independent_path = if_else(is.na(independent_path), "Ind", independent_path)) %>% 
    left_join(vape_res, by = c("leadingEdge" = "ENTREZID"))
  
  print(head(long_tab))
  long_tab <- long_tab %>% 
    dplyr::select(pathway, independent_path, leadingEdge, 
                  gene_name, log2FoldChange, regulation, pvalue, padj.y, padj.x)
  
  return(long_tab)
}
#Get supplementary tables
kegg_res_full <- sup_table_maker(kegg_ind)

react_res_full <- sup_table_maker(react_ind)

go_res_full <- sup_table_maker(go_ind)

colnames(kegg_res_full) <- c("Pathway", "Independent Pathway", "leadingEdge", "gene_name", "log2FoldChange", "regulation", "gene_level_pval", "gene_level_fdr", "path_level_fdr")
colnames(react_res_full) <- c("Pathway", "Independent Pathway", "leadingEdge", "gene_name", "log2FoldChange", "regulation", "gene_level_pval", "gene_level_fdr", "path_level_fdr")
colnames(go_res_full) <- c("Pathway", "Independent Pathway", "leadingEdge", "gene_name", "log2FoldChange", "regulation", "gene_level_pval", "gene_level_fdr", "path_level_fdr")

#Write out kegg results
list_of_datasets <- list("pathway_level" = kegg_res_supplement, "gene_level" = kegg_res_full)
writexl::write_xlsx(list_of_datasets, here("DataProcessed/gsea/kegg_full_res_2022_10_20.xlsx"))

#Write out react results
list_of_datasets <- list("pathway_level" = react_res_supplement, "gene_level" = react_res_full)
writexl::write_xlsx(list_of_datasets, here("DataProcessed/gsea/react_full_res_2022_10_20.xlsx"))

#Write out go results
list_of_datasets <- list("pathway_level" = go_res_suplement, "gene_level" = go_res_full)
writexl::write_xlsx(list_of_datasets, here("DataProcessed/gsea/go_full_res_2022_10_20.xlsx"))

write_csv(go_res_full, here("DataProcessed/gsea/go_gene_level_res_2023_01_19.csv"))
