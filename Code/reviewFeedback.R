

rm(list = ls())

library(tidyverse)
library(ENmix)

`%notin%` <- Negate(`%in%`)

setwd("P:/BRANCHES/Pulmonary/Sharma_Vaping/")

###################
# 1. We could do a correlation analysis between cg02123174 & REXO1 only.
###################

#Read in the m-values
bvals <- read_tsv("DataProcessed/methylation/QC/methylation_betas_final_2022_09_27.txt") %>% column_to_rownames("CpG_Site") 
mvals <- B2M(bvals)
vstcnts  <- as.data.frame(read_csv("DataProcessed/rna_seq/ruv/RUV_k2_norm_counts_2022_10_13.csv"))
metadata <- as.data.frame(read_csv("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv"))

mvals01 <- t(mvals["cg02123174",])
meta01 <- merge(mvals01, metadata, by.x = 0, by.y = "sentrix_name")
genes01 <- t(vstcnts[grepl("ENSG00000079313", vstcnts$gene),])
meta02 <- merge(meta01, genes01, by.x = "rna_id", by.y = 0)

meta01[which(meta01$sid %notin% meta02$sid),c("sid", "rna_id", "cg02123174")]
genes01[which(rownames(genes01) %notin% meta01$rna_id),]

summary(lm(`1277` ~ cg02123174 + age + recruitment_center + sex_lab + vape_6mo_lab, data = meta02))

plot(meta02$cg02123174, as.numeric(meta02$`1277`), xlab = "M-Values of cg02123174", ylab = "Normalized REXO1 Counts", main = "Pearson's Correlation = -0.06\np-value = 0.7")

###################
# 2. For the DMR results, we could do something similar with the DMRs with their closest gene – but would need to figure out some details.
###################

dmr.res1 <- read_csv("DataProcessed/methylation/DMR/results/combp_format_res.csv")
dmr.res2 <- read_csv("DataProcessed/methylation/DMR/results/dmr_res_ENmix.csv")
dmr.res3 <- read_csv("DataProcessed/methylation/DMR/results/dmr_res_python.csv")

dmr.res1 %>% filter(sidak < 0.1) %>% arrange(chrom)
dmr.res2 %>% filter(sidak < 0.1) %>% arrange(chr)
dmr.res3 %>% filter(z_sidak_p < 0.1) #%>% arrange(`#chrom`)

#dmRegions <- dmr.res %>% filter(sidak < 0.1 & is.na(gene_name) == F) 

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(org.Hs.eg.db)
library(stringr)
library(biomaRt)

## Read in the meth results
meth_res <- read_csv("DataProcessed/methylation/EWAS/results/full_res_bacon_2022_12_01.csv")

## Access annotation data
data("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
annotation.table <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

## Convert to Data Frame for easier Access
mapping <- data.frame("Name" = annotation.table$Name,
                      "RefSeq" = annotation.table$UCSC_RefGene_Name,
                      "chrom" = annotation.table$chr,
                      "pos" = annotation.table$pos, 
                      "strand" = annotation.table$strand,
                      "island_name" = annotation.table$Islands_Name,
                      "relation_to_island" = annotation.table$Relation_to_Island)

## Merge With EWAS results 
cpg.island.map <- left_join(meth_res, mapping, by = c("CpG_Site" = "Name"))

## Format fields to match bed file 
cpg.island.map <- cpg.island.map %>% 
  mutate(start = pos,
         end = start + 51,
         chrom = gsub("chr", "", chrom),
         chr = factor(chrom, levels = paste0(c(seq(1:22), "X", "Y"))))

#Read back in the results from comb-p
dmr.res <- read_tsv("DataProcessed/methylation/DMR/results/dmr_res_py.bed", show_col_types = F) %>% 
  janitor::clean_names() %>% 
  dplyr::rename(chrom = number_chrom)
#Produce key file for dmp res
dmp.key <- cpg.island.map %>% 
  dplyr::select(chrom, pos, CpG_Site, Estimate, pval.bacon) %>%
  drop_na(chrom) %>% 
  mutate(chrom = paste0("chr", chrom))

#Find which cpgs that fall within each DMR and sum up if they are pos or neg methylated
get_meth_dir <- function(dmr, dmp) {
  tmp_res <- matrix(ncol = 3, nrow = dim(dmr)[1])
  colnames(tmp_res) <- c("cpg_sites","n_meth_pos", "n_meth_neg")
  
  for (i in 1:dim(dmr)[1]){
    chr <- dmr$chrom[[i]]
    start <- dmr$start[[i]]
    end <- dmr$end[[i]]
    minp <- dmr$min_p[[i]]
    
    temp <- dmp[dmp$chrom == chr & between(dmp$pos, start, end),]
    cpg_list <- paste(temp$CpG_Site, collapse = ";")
    tmp_res[i,] <- c(cpg_list, sum(temp$Estimate > 0), sum(temp$Estimate < 0))
  }
  return(cbind(dmr, tmp_res))
}

#Running the function 
dmr.full.anno <- get_meth_dir(dmr.res, dmp.key) %>% 
  mutate(n_meth_pos = as.numeric(n_meth_pos),
         n_meth_neg = as.numeric(n_meth_neg)) %>% 
  arrange(z_sidak_p)

dmr.full.anno.filter <- dmr.full.anno %>% filter(z_sidak_p < 0.1) %>% arrange(z_sidak_p) %>% as.data.frame() %>% 
  dplyr::select(chrom, start, end, min_p, n_probes, z_p, z_sidak_p, ref_gene_name, cpg_island_ext_feature, cpg_sites)


i <- 1
j <- 1
k <- 1

# doesn't have: EIPR1, PPP1R3A, SLITRK5, NOC4L
vstcnts %>% filter(gene == "")

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")

key <- getBM(attributes = c('ensembl_gene_id_version', 'external_gene_name'),
             filters = 'ensembl_gene_id_version',
             values = vstcnts$gene, 
             mart = ensembl)
cnts <- merge(key, vstcnts, by.x = "ensembl_gene_id_version", by.y = "gene")
genes <- dmr.full.anno.filter$ref_gene_name
genes <- unique(unlist(str_split(genes, ";")))
cnts <- cnts %>% filter(cnts$external_gene_name  %in% genes)
rownames(cnts) <- cnts$external_gene_name
cnts

# doesn't contain: PPP1R3A, SLITRK5, NBL1, MICOS10-NBL1, HDAC4

dmr.full.anno.filter2 <- dmr.full.anno.filter[grepl(paste0(cnts$external_gene_name, collapse = "|"), dmr.full.anno.filter$ref_gene_name),]

# only have information for one genes
dmr.full.anno.filter2[grepl("MICOS10-NBL1", dmr.full.anno.filter2$ref_gene_name), "ref_gene_name"] <- "NBL1"
dmr.full.anno.filter2[grepl("LOC100288123", dmr.full.anno.filter2$ref_gene_name), "ref_gene_name"] <- "REXO1"

#just using to figure out the number of rows i need
s <- 0

outDf <- data.frame("cpg"  = rep(NA, 75), 
                    "gene" = rep(NA, 75), 
                    "estimate"  = rep(NA, 75), 
                    "p"    = rep(NA, 75))

for(i in 1:nrow(dmr.full.anno.filter2)){
  genes <- dmr.full.anno.filter2$ref_gene_name[i]
  cpgs <- dmr.full.anno.filter2$cpg_sites[i]
  
  if(grepl(";", genes)){genes <- str_split(genes, ";", simplify = T)}
  if(grepl(";", cpgs)){cpgs <- str_split(cpgs, ";", simplify = T)}
  
  for(j in 1:length(cpgs)){
    for(k in 1:length(genes)){
      mvals01 <- t(mvals[cpgs[j],])
      meta01 <- merge(mvals01, metadata, by.x = 0, by.y = "sentrix_name")
      genes01 <- t(cnts[grepl(genes[k], cnts$external_gene_name),])
      if(ncol(genes01) == 0){
        print(i)
        print(j)
        print(k)
      }
      meta02 <- merge(meta01, genes01, by.x = "rna_id", by.y = 0)
      
      out <- summary(lm(get(genes[k]) ~ get(cpgs[j]) + age + recruitment_center + sex_lab + vape_6mo_lab, data = meta02))
    
      s <- s + 1 
    
      outDf$cpg[s]  <- cpgs[j]
      outDf$gene[s] <- genes[k]
      outDf$estimate[s]  <- out$coefficients["get(cpgs[j])", "Estimate"]
      outDf$p[s]    <- out$coefficients["get(cpgs[j])", "Pr(>|t|)"]
      
      if(out$coefficients["get(cpgs[j])", "Pr(>|t|)"] < 0.1){
        png(filename = paste0("P:/BRANCHES/Pulmonary/Sharma_Vaping/figures/ReviewerFigures/lmPlot_", cpgs[j], "_", genes[k], ".png"))
        plot(meta02[,cpgs[j]], as.numeric(meta02[,genes[k]]), 
             xlab = paste0("M-Values of ", cpgs[j]), 
             ylab = paste0("Normalized ", genes[k], " Counts"), 
             main = paste0("Linear Model Estimate = ", round(out$coefficients["get(cpgs[j])", "Estimate"], digits = 0), "\np-value = ",  round(out$coefficients["get(cpgs[j])", "Pr(>|t|)"], 2)))
        dev.off()
      }
    }
  }
}

# remove the duplicates
outDf <- outDf %>% distinct(cpg, gene, .keep_all = T)

# calculate FDR
outDf$FDR <- p.adjust(outDf$p, method = "BH")

outDf %>% arrange(p) %>% filter(p < 0.05)
outDf %>% arrange(p) %>% filter(p < 0.1)

write.csv(outDf %>% arrange(p), row.names = F, 
          file = "P:/BRANCHES/Pulmonary/Sharma_Vaping/DataProcessed/reviewerAdditions/lmResults_for_DMR.csv")

