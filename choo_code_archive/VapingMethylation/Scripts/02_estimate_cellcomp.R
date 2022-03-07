cd /beevol/home/borengas/analyses/methylation


module load R/3.6.1
R


# reload targets 01, but using minfi infrastructure for compatbility
# load idats via minfi
rgset_minfi <-
  read.metharray.exp(
    targets = targets %>%
      transmute(Basename = file_path_prefix))


# experimental hub: EPIC blood dataset
library(ExperimentHub)  
hub <- ExperimentHub()  
query(hub, "FlowSorted.Blood.EPIC")  

FlowSorted.Blood.EPIC <- hub[["EH1136"]]  


library(FlowSorted.Blood.EPIC)
# run & save cell counts
cellcomp <-
  estimateCellCounts2(rgset_minfi,
                    compositeCellType = "Blood",   
                    processMethod = "preprocessNoob",  
                    cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),  
                    referencePlatform = "IlluminaHumanMethylationEPIC",  
                    IDOLOptimizedCpGs = IDOLOptimizedCpGs,
                    lessThanOne = T)


cellcomp$counts %>%
  as_tibble(rownames = "Plate_ID") %>%
  write_tsv(.,
            path = paste0("./Output/", date_export, "_MethylationDat/cellcomp.tsv"))