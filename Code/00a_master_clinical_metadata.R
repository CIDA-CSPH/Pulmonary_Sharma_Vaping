####################Load Libraries########################
library(tidyverse)
library(here)
library(readxl)
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
id_relate_meth <- read_tsv(file = here("DataRaw/subject_ids/20210401_methylation_coreID_to_PID.txt"), col_names = T) %>% clean_names()
id_relate_rna <- read_tsv(file = here("DataRaw/subject_ids/20201216_coreID_to_PID.txt"), col_names = T) %>% clean_names()
metadata_unjoined <- read_csv(file = here("DataProcessed/clinical_metadata/table1_clean_data_2022_08_22.csv"))

#Join Sample ID to Subject ID
metadata_joined <- id_relate_meth %>% 
  mutate(methylation_id = str_pad(new_id, 2, "left", "0") %>% paste0("Sample", .)) %>% 
  select(-new_id) %>% 
  full_join(metadata_unjoined, by = "sid")

metadata_joined <- id_relate_rna %>% 
  mutate(rna_id = str_pad(new_id, 2, "left", "0") %>% paste0("Sample", .)) %>% 
  select(-new_id) %>% 
  full_join(metadata_joined, by = "sid")

Sentrix_ID_to_Meta_ID <- full_join(metadata_joined, samplesheet %>% select(sample_join, sentrix_name, `Sentrix Barcode`, `Sentrix Position`), 
                                   by = c("methylation_id" = "sample_join"))

write_csv(Sentrix_ID_to_Meta_ID, here("/Users/hawkjona/Documents/CIDA/Sharma_Vaping_Local.nosync/DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv"))
