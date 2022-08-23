####################Load Libraries########################
library(tidyverse)
library(here)
library(sesame) #ver 1.15.6
library(readxl)
library(parallel)
library(randomForest)
library(janitor)

################## Setup #################
#Read in and join methylation metadata
samplesheet <- left_join(
  read_csv(here("DataRaw/methylation/Metadata/SSharma_48_samplesheet_05202021.csv"), skip = 7) %>%
    select(-`Sample_Group`, -`Pool_ID`),
  read_excel(here("DataRaw/methylation/Metadata/Sharma_48_5192021Controls.xlsx")) %>%
    mutate(`Sample Name` = as.numeric(`Sample Name`)),
  by = c("Sample_Name" = "Sample Name"))

samplesheet <- samplesheet %>% 
  mutate(sample_join = paste0("Sample", str_pad(as.character(Sample_Name), 2, pad = "0")),
         sentrix_name = paste0(`Sentrix Barcode`, "_", `Sentrix Position`))

#Read in and join Clinical Metadata
#Load 
id_relate <- read_tsv(file = here("DataRaw/methylation/Metadata/20210401_methylation_coreID_to_PID.txt"), col_names = T) %>% clean_names()
metadata_unjoined <- read_csv(file = here("DataProcessed/metadata/table1_clean_data_2022_08_22.csv"))

#Join Sample ID to Subject ID
metadata_joined <- id_relate %>% 
  mutate(new_id = str_pad(new_id, 2, "left", "0") %>% paste0("Sample", .)) %>% 
  left_join(metadata_unjoined, by = "sid")

Sentrix_ID_to_Meta_ID <- left_join(metadata_joined, samplesheet %>% select(sample_join, sentrix_name), 
                                   by = c("new_id" = "sample_join"))

################## Read in IDAT data #################
folder_raw_dat <- here("DataRaw/methylation/RawIdat/")

methylation_raw <- lapply(searchIDATprefixes(folder_raw_dat),
                          readIDATpair)

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
metadata_pred_sex <- left_join(Sentrix_ID_to_Meta_ID, pred_sex_df, by = "sentrix_name")

metadata_all_sex <- left_join(metadata_pred_sex, chr_intensity, by = "sentrix_name")

metadata_all_sex[,18:46] <- lapply(metadata_all_sex[,18:46], as.numeric)

# write_csv(metadata_all_sex, here("DataProcessed/methylation/metadata_all_sex_2022_08_22.csv"))


