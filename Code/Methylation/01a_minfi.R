####################Load Libraries########################
library(tidyverse)
library(here)
library(kableExtra)
library(minfi)
library(readxl)
library(parallel)
library(randomForest)
library(janitor)


#################### Set Base Directory ########################

baseDir <- here("DataRaw/methylation/RawIdat")

list.files(baseDir)

#################### Read in Sample Data ########################
################## Setup #################
#Read in and join metadata
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

#################### Read in Sample Data ########################

RGset <- read.metharray.exp(targets = targets)

mdsPlot(RGset, sampNames = RGset$Sample_Name)
