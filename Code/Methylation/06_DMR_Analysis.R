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

data("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
annotation.table <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

mapping <- data.frame("Name" = annotation.table$Name,
                      "RefSeq" = annotation.table$UCSC_RefGene_Name,
                      "chrom" = annotation.table$chr,
                      "pos" = annotation.table$pos, 
                      "strand" = annotation.table$strand,
                      "island_name" = annotation.table$Islands_Name,
                      "relation_to_island" = annotation.table$Relation_to_Island)

dim(mapping)


#Unlist genes and islands
genelist <- strsplit(mapping$RefSeq, split = ";")
names(genelist) <- mapping$island_name
islands <- unlist(genelist)
islands <- islands[names(islands) != ""]



#Create Key file for islands to genes
islands.genes.map <- data.frame(island_name = names(islands),
                                gene_name = islands) %>% 
  distinct()

## Make key file for CpG's to islands
cpg.island.map <- mapping %>% 
  mutate(island_name = gsub(".*:", "", island_name),
         start = sapply(strsplit(island_name, "-"), `[`, 1),
         end = sapply(strsplit(island_name, "-"), `[`, 2))



## Join methylation results
CpG.Mapped <- left_join(meth_res, cpg.island.map, by = c("CpG_Site" = "cpg")) %>% 
  drop_na(pos) %>% 
  dplyr::select(-c(p.value, fdr, island_name))

## Comb-p formatting
dmr.full <- CpG.Mapped %>% 
  mutate(chr = gsub("chr", "", chr)) %>% 
  dplyr::select(chr, pos, end, pval.bacon, CpG_Site) %>% 
  dplyr::rename(p = pval.bacon,
                probe = CpG_Site,
                start = pos)

#Make Chromosome a factor
dmr.full$chr <- factor(dmr.full$chr, levels = paste0(c(seq(1:22), "X", "Y")))

#Remove leading 0 from probe names
clean.cpg <- str_remove(dmr.full$probe, "cg") %>% str_remove("^0+")
dmr.full$probe <- paste0("cg", clean.cpg)

# Convert start and end to numberic
dmr.full$start <- as.integer(dmr.full$start)
dmr.full$end <- as.integer(dmr.full$end)

## Run Comb-p in wrapper
combp(dmr.full)

## Read back in the results
dmr.res <- read_csv(here("DataProcessed/methylation/results/DMR/resu_combp.csv"))

dmr.res$name <- paste0("chr",dmr.res$chr, ":", dmr.res$start, "-", dmr.res$end)

dmr.full %>% 
  filter(start %in% dmr.res$start)

dmr.res %>% 
  filter(name %in% islands.genes.map$island_name)


dat=simubed()
names(dat)
#seed=0.1 is only for demonstration purpose, it should be smaller than 0.05 or 0.01 in actual study.
combp(data=dat,seed=0.1)
