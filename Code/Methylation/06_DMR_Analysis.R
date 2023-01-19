## ---------------------------
## Script name: DMR Analysis
##
## Purpose of script: Complete full DMR analysis
##
## Author: Trent Hawkins
##
## Date Created: 2023-01-18
## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(here)
require(ENmix)
require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
require(org.Hs.eg.db)
require(missMethyl)

## Read in the meth results
meth_res <- read_csv(here("DataProcessed/methylation/results/full_res_bacon_2022_12_01.csv")) 



## Filter Significant results
meth_res.sig <- meth_res %>% filter(fdr.bacon < 0.05)


## Get Illumina Methylation
anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

## Get needed annotations
anno_df <- data.frame(cpg = anno@listData$Name,
                      chr = anno@listData$chr,
                 island_name = anno@listData$Islands_Name)

## Make key to merge with results
df <- anno_df %>% 
  mutate(island_name = gsub(".*:", "", island_name),
         start = sapply(strsplit(island_name, "-"), `[`, 1),
         end = sapply(strsplit(island_name, "-"), `[`, 2))

CpG.Mapped <- left_join(meth_res, df, by = c("CpG_Site" = "cpg")) %>% 
  drop_na(start) %>% 
  dplyr::select(-c(p.value, fdr, island_name))

dmr.full <- CpG.Mapped %>% 
  mutate(chr = gsub("chr", "", chr)) %>% 
  dplyr::select(chr, start, end, pval.bacon, CpG_Site) %>% 
  dplyr::rename(p = pval.bacon,
                probe = CpG_Site)

dmr.test <- dmr.full %>% 
  filter(chr == 1)

combp(dmr.full, nCores = detectCores())
