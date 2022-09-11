####################Load Libraries########################
library(tidyverse)
library(here)
library(kableExtra)
library(sesame) #ver 1.15.6
library(readxl)
library(parallel)
library(randomForest)
library(janitor)
library(data.table)

kablize <- function(tab, digits = 3) {
  kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}

format_num <- function(number, digits = 0, format = "f") {
  formatC(number, digits = digits, format = format, big.mark = ",")
}
################## Setup #################

#Read in master clinical metadata
clin_metadata <- read_csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv"))

#Read in and join methyl metadata
samplesheet <- left_join(
  read_csv(here("DataRaw/methylation/Metadata/SSharma_48_samplesheet_05202021.csv"), skip = 7) %>%
    select(-`Sample_Group`, -`Pool_ID`),
  read_excel(here("DataRaw/methylation/Metadata/Sharma_48_5192021Controls.xlsx")) %>%
    mutate(`Sample Name` = as.numeric(`Sample Name`)),
  by = c("Sample_Name" = "Sample Name"))

samplesheet <- samplesheet %>% 
  mutate(sample_join = paste0("Sample", str_pad(as.character(Sample_Name), 2, pad = "0")),
         sentrix_name = paste0(`Sentrix Barcode`, "_", `Sentrix Position`))

# Getting Filename prefixes for readIDATpairs
folder_raw_dat <- here("DataRaw/methylation/RawIdat/")

targets <- tibble(file_path =
                    list.files(path = folder_raw_dat, pattern = "*.idat|*.IDAT",
                               full.names = F, recursive = T)) %>%
  mutate(`Sentrix Barcode` = str_split(file_path, pattern = "/|_") %>% map_chr( ~ .[4]), # chip
         `Sentrix Position` = str_split(file_path, pattern = "/|_") %>% map_chr( ~ .[5])) %>% # row
  left_join(samplesheet,
            by = c("Sentrix Barcode", "Sentrix Position")) %>%
  mutate(file_path_prefix = gsub("_Red.idat|_Grn.idat", "", file_path) %>%
           paste0(folder_raw_dat, .)) %>%
  filter(!duplicated(file_path_prefix))

targets <- targets %>% 
  dplyr::rename(Basename = file_path_prefix)

############## SeSAMe Pipeline By Hand (FINICKY B/C Parallelization!!!!)########
# ################## Read in IDAT data #################
# methylation_raw <- lapply(searchIDATprefixes(folder_raw_dat),
#                                           readIDATpair)
# 
# ################## Run openSesame pipeline step-by-step #################
# #qualityMask
# methylation_mask <- BiocParallel::bplapply(methylation_raw, qualityMask)
# #Infer Infinium I Channel
# methylation_infer <- BiocParallel::bplapply(methylation_mask, inferInfiniumIChannel)
# #Dye-Bias Correction
# methylation_db_corr <- BiocParallel::bplapply(methylation_infer, dyeBiasNL)
# #detection p-value masking
# methylation_detect_mask <- BiocParallel::bplapply(methylation_db_corr, pOOBAH)
# #Background subtraction
# methylation_background <- BiocParallel::bplapply(methylation_detect_mask, noob)
# #getBetas
# betas <- do.call(cbind, BiocParallel::bplapply(methylation_background, getBetas)) %>% as.data.frame()

################ get Betas with openSesame Pipeline #########
#Noob only
#betas <- openSesame(folder_raw_dat) %>% as.data.frame()

#Noob + BMIQ
betas <- openSesame(folder_raw_dat, prep = "QCDPBM") %>% as.data.frame()
################## Drop probes with NA for > 10% of Samples #################

cpg_drop_na <- function(betas) {
  # how many NAs for each CPG Site
  cpg_na_count <-  apply(betas, 1, function(x) {sum(is.na(x))})
  
  #Get the CPG sites that have NA's for > 10% of samples
  cpg_drop <- cpg_na_count > 0.1*ncol(betas)
  
  #Get Betas with the removed CPGs
  betas_dropped <- betas[!cpg_drop,]
  
  return(betas_dropped)
}

betas_drop <- cpg_drop_na(betas)

exp_dropped <- dim(betas)[1] - dim(betas_drop)[1] - 105454 #105454 is number masked by qualityMask()

#Turn CPG into it's own column
betas <- betas %>%
  mutate(CpG_Site = rownames(.)) %>%
  select(CpG_Site, everything())

betas_drop <- betas_drop %>%
  mutate(CpG_Site = rownames(.)) %>%
  select(CpG_Site, everything())

################# Write out Betas #####################

#write_tsv(betas_drop, file = here("DataProcessed/methylation/methylation_betas_BMIQ.txt"))

#write_tsv(betas, file = here("DataProcessed/methylation/methylation_betas_full_BMIQ.txt"))


############## Convert to M Values #####################

mvals <- betas %>% 
  select(-CpG_Site) %>% 
  apply(., 2, BetaValueToMValue) %>% 
  as.data.frame() %>% 
  mutate(CpG_Site = rownames(.)) %>% 
  select(CpG_Site, everything())

mvals_drop <- cpg_drop_na(mvals)

################# Write out M-Values #####################
 #write_tsv(mvals_drop, file = here("DataProcessed/methylation/methylation_mvals_BMIQ.txt"))
# 
 #write_tsv(mvals, file = here("DataProcessed/methylation/methylation_mvals_full_BMIQ.txt"))

################ Get Betas using only Autosomes #####################
betas <- read_tsv(here("DataProcessed/methylation/methylation_betas_full_BMIQ.txt"))

get_auto_vals <- function(betas, get_value = "Betas"){
  #Load in raw idats
  methylation_raw <- lapply(searchIDATprefixes(folder_raw_dat), readIDATpair)
  
  #Get annotation
  autosomes <- sesameData_getAutosomeProbes("EPIC")
  
  cpg_sites <- names(autosomes)
  
  probe_annos <- autosomes %>% as.data.frame(row.names = names(autosomes))
  
  autosomal_probes <- probe_annos[!probe_annos$seqnames %in% c("chrX", "chrY", "chrM"),]
  
  #Filter for only autosomal probes
  methyl_raw_auto <- lapply(methylation_raw, function(x) {
    x %>% 
      filter(Probe_ID %in% rownames(autosomal_probes))
  })
  #Noob Only 
  betas_auto <- openSesame(methyl_raw_auto, prep = "QCDPBM") %>% as.data.frame()
  return(betas_auto)
  }

betas_auto <- get_auto_vals(betas)

#Turn CPG into it's own column
betas_auto <- betas_auto %>%
  mutate(CpG_Site = rownames(.)) %>%
  select(CpG_Site, everything())

betas_auto_drop <- cpg_drop_na(betas_auto)

############### Write out Betas #######################

#With NA's Dropped
write_tsv(betas_auto_drop, here('DataProcessed/methylation/autosomal_betas_BMIQ.txt'))

#Full DF
write_tsv(betas_auto, here('DataProcessed/methylation/autosomal_betas_full_BMIQ.txt'))



############### Convert to Mvals #######################
#betas_auto <- read_tsv(here('DataProcessed/methylation/autosomal_betas.txt'))

m_auto <- betas_auto %>%
  select(-CpG_Site) %>% 
  as.matrix() %>% 
  BetaValueToMValue() %>%
  as.data.frame() %>%
  mutate(CpG_Site = betas_auto$CpG_Site) %>% 
  select(CpG_Site, everything())

m_auto_drop <- cpg_drop_na(m_auto)
############### Write out mvals #######################
# #With NA's Dropped
# write_tsv(m_auto_drop, here('DataProcessed/methylation/autosomal_mvals.txt'))
# 
# #Full DF
# write_tsv(m_auto, here('DataProcessed/methylation/autosomal_mvals_full.txt'))

