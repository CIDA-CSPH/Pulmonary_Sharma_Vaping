## ---------------------------
## Script name: Map CpG to Genome
##
## Purpose of script: Targeted Methylation Analysis
##
## Author: Trent Hawkins
##
## Date Created: 2023-01-04
## ---------------------------
## set working directory for Mac and PC

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(here)
require(org.Hs.eg.db)
require(minfi)
require(readxl)
require(missMethyl)
require(clusterProfiler)
require(openxlsx)


# Needed function from missMethyl source code -----------------------------

.getFlatAnnotation <- function(array.type=c("450K","EPIC"),anno=NULL)
  # flatten 450k or EPIC array annotation
  # Jovana Maksimovic
  # 18 September 2018
  # Updated 18 September 2018
  # Modified version of Belida Phipson's .flattenAnn code
{
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  
  # get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$Name),]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  
  flat$cpg<- substr(rownames(flat),1,10)
  
  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma::alias2SymbolTable(flat$symbol))
  
  #eg <- toTable(org.Hs.egSYMBOL2EG)
  eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys=AnnotationDbi::keys(org.Hs.eg.db), 
                                               columns=c("ENTREZID","SYMBOL"), 
                                               keytype="ENTREZID"))
  colnames(eg) <- c("gene_id","symbol")
  
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  return(flat)
}

# Read in RNASeq Results, Map to ENTREZ, and extract sig genes --------------------------------------------------

## Read in results
RNAseq.res <- read_csv(here("DataProcessed/rna_seq/differential_expression/full_analysis/full_vape_res_2022_10_13.csv"))

## Get significant genes
RNAseq.sig <- RNAseq.res %>% 
  filter(padj < 0.05)


# Read Methylation Results and extract CpG Sites mapped to sig Genes --------

## Read in the meth results
meth_res <- read_csv(here("DataProcessed/methylation/EWAS/results/full_res_bacon_2022_12_01.csv"))

## Genes to CpG Key
genes_to_cpg <- .getFlatAnnotation(array.type = "EPIC")

## Filter for CpGs associated with significant genes
sigGene.map <- genes_to_cpg %>% 
  dplyr::filter(symbol %in% RNAseq.sig$gene_name)

## Get unique CpG loci
sigGene.cpgs <- unique(sigGene.map$cpg)

## Subset sig loci from meth results
methRes.target <- meth_res %>% 
  dplyr::filter(CpG_Site %in% sigGene.cpgs)

## Re-adjust the p-values
methRes.target$fdr.target <- p.adjust(methRes.target$pval.bacon, method = "fdr")

## Subset new sig cpgs
sig.target.cpgs <- methRes.target %>% 
  dplyr::filter(fdr.target < 0.1)

## Map CpGs back to genes
sigGene.map %>% 
  dplyr::filter(cpg %in% sig.target.cpgs$CpG_Site)

## Create dataframe for final results
sigCites.final <- left_join(sigGene.map, methRes.target %>% dplyr::select(CpG_Site, Estimate, `Std. Error`, `t value`, pval.bacon, fdr.target), by = c("cpg" = "CpG_Site")) %>% 
  drop_na(fdr.target) %>% 
  dplyr::select(cpg, symbol, Estimate, `Std. Error`, `t value`, pval.bacon, fdr.target) %>% 
  distinct()

## Only results we want to report to the investigator
sigCites.report <- sigCites.final %>%
  left_join(., RNAseq.res, by = c("symbol" = "gene_name")) %>% 
  dplyr::mutate(gene.reg = if_else(log2FoldChange > 0, '+', '-')) %>% 
  filter(fdr.target < 0.2) %>% 
  dplyr::arrange(fdr.target) %>% 
  distinct() %>% 
  dplyr::select(cpg, symbol, Estimate, `Std. Error`, `t value`, pval.bacon, fdr.target, gene.reg) %>% 
  dplyr::rename("CpG Site" = cpg,
                "Symbol" = symbol,
                "Estimate (M)" = Estimate,
                "pval" = pval.bacon,
                "FDR" = fdr.target)

## Get records for related genes
sigCites.generelate <- RNAseq.res %>% 
  dplyr::filter(gene_name %in% sigCites.report$Symbol) %>% 
  dplyr::rename("EnsemblID" = ensg,
                "Symbol" = gene_name,
                "Std. Error" = lfcSE,
                "t-value" = stat,
                "pval" = pvalue,
                "FDR" = padj) %>% 
  dplyr::select(-c(gene_type, baseMean))

## Write out the results
#Write out go results
list_of_datasets <- list("DMP" = sigCites.report, "Associated Genes" = sigCites.generelate)
writexl::write_xlsx(list_of_datasets, here("DataProcessed/methylation/Targeted/results/targetedDMP_results_reporting_2023_01_26.xlsx"))
