####################Load Libraries########################
library(tidyverse)
library(here)
library(sesame) #ver 1.15.6
library(readxl)
library(parallel)
library(randomForest)
library(janitor)

################## Setup #################
#Read in master clinical metadata
clin_metadata <- read_csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv"))

################## Read in IDAT data #################
folder_raw_dat <- here("DataRaw/methylation/RawIdat/")

#Read in to unprocessed SigDF
methylation_raw <- lapply(searchIDATprefixes(folder_raw_dat), readIDATpair)

############## Get Sex Info and Compare #####################
#Get predicted Sex
predicted_sex <- lapply(methylation_raw, inferSex)

pred_sex_df <- stack(predicted_sex)

colnames(pred_sex_df) <- c("pred_sex", "sentrix_name")

#Get Chromosome intensities
chr_intensity <- lapply(methylation_raw, getSexInfo) 

chr_intensity <- data.frame(t(data.frame(chr_intensity))) %>%
  dplyr::mutate("sentrix_name" = rownames(.) %>% gsub("X", "", .))


#Join with metadata
metadata_pred_sex <- left_join(clin_metadata, pred_sex_df, by = "sentrix_name")

metadata_all_sex <- left_join(metadata_pred_sex, chr_intensity, by = "sentrix_name")

#write_csv(metadata_all_sex, here("DataProcessed/methylation/metadata_all_sex_yyyy_mm_dd.csv"))


