####################Load Libraries########################
library(tidyverse)
library(here)
library(minfi)
library(readxl)
library(bumphunter)
library(RColorBrewer)

#################### Set Base Directory ########################

baseDir <- here("DataRaw/methylation/RawIdat")

list.files(baseDir)

#################### Read in Sample Data ########################
################## Setup #################
#Read in clinical metadata
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

#join metadata back to targets
targets <- left_join(targets, clin_metadata %>% 
                       select(methylation_id,sid,sex_lab, age, recruitment_center, vape_6mo_lab),
                     by = c("sample_join" = "methylation_id"))


#################### Read in Sample Data ########################

RGset <- read.metharray.exp(targets = targets)

getManifest(RGset)

minfi::getProbeInfo(RGset)

#################### Plotting and QC ########################
#Beta Density Plots
densityPlot(RGset, sampGroups = targets$sex_lab, main = "Beta", xlab = "Beta")
densityPlot(RGset, sampGroups = targets$vape_6mo_lab, main = "Beta", xlab = "Beta")
densityBeanPlot(RGset, sampGroups = targets$sex_lab, sampNames = targets$sid)
densityBeanPlot(RGset, sampGroups = targets$vape_6mo_lab, sampNames = targets$sid)

##Control Plots
#Negative Controls
controlStripPlot(RGset, controls = c("NEGATIVE"), sampNames = targets$Sample_Name)
#Bisulfite Conversion Controls
controlStripPlot(RGset, controls = c("BISULFITE CONVERSION I"), sampNames = targets$Sample_Name)

################ Detection P-Values ##########################
detP <- detectionP(RGset)
failed <- detP > 0.05
numfail.col = colMeans(failed)
which.max(numfail.col)
badProbes <- rowMeans(detP) >= 0.05
sum(badProbes)

################# Get Annotation ############################

################ Normalization and MDS Plots (Whole Genome) ##################
mset = preprocessRaw(RGset)
set.seed(404)
msetNoob <- preprocessNoob(RGset)

annotation <- getAnnotation(msetNoob)
autosomes <- annotation[!annotation$chr %in% c("chrX", "chrY"),]
msetNoob_auto <- msetNoob[featureNames(msetNoob) %in% row.names(autosomes),]
beta_noob_auto <- getBeta(msetNoob_auto)
m_noob_auto    <- getM(msetNoob_auto)
colnames(m_noob_auto) <- sampleNames(msetNoob)
totalProbes2 <- prettyNum(dim(beta_noob_auto)[1], big.mark = ",")



# num_1000 <- 1000
# num_10000 <- 10000
# sex_pal <- c("deepskyblue", "deeppink3")
# center_colors <- c("cyan3", "chartreuse3", "darkorange")
# 
# #Raw data
# par(mfrow = c(2,3))
# mdsPlot(mset, sampGroups = targets$sex_lab, sampNames = targets$Sample_Name,
#         numPositions=num_1000, ylim = c(-10,10), main = "MDS - Raw data - Sex -
#  1000 sites", pal = sex_pal)
# 
# mdsPlot(mset, sampGroups = targets$vape_6mo_lab, sampNames = targets$Sample_Name,
#         numPositions=num_1000, ylim = c(-10,10), main = "MDS - Raw data - Vape Status -
#  1000 sites")
# 
# mdsPlot(mset, sampGroups = targets$recruitment_center, sampNames = targets$Sample_Name,
#         numPositions=num_1000, ylim = c(-10,10), pal = center_colors, main = "MDS - Raw data - Recruitment Center -
#  1000 sites")
# 
# mdsPlot(mset, sampGroups = targets$sex_lab, sampNames = targets$Sample_Name,
#         numPositions=num_10000, ylim = c(-10,10), main = "MDS - Raw Data - Sex -
#  10000 sites", pal = sex_pal)
# 
# mdsPlot(mset, sampGroups = targets$vape_6mo_lab, sampNames = targets$Sample_Name,
#         numPositions=num_10000, ylim = c(-10,10), main = "MDS - Raw Data - Vape Status -
#  10000 sites")
# 
# mdsPlot(mset, sampGroups = targets$recruitment_center, sampNames = targets$Sample_Name,
#         numPositions=num_10000, ylim = c(-10,10), pal = center_colors, main = "MDS - Raw data - Recruitment Center -
#  10000 sites")
# 
# #Noob data
# par(mfrow = c(2,3))
# mdsPlot(msetNoob, sampGroups = targets$sex_lab, sampNames = targets$Sample_Name,
#         numPositions=num_1000, ylim = c(-10,10), main = "MDS - Noob - Sex -
#  1000 sites", pal = sex_pal)
# 
# mdsPlot(msetNoob, sampGroups = targets$vape_6mo_lab, sampNames = targets$Sample_Name,
#         numPositions=num_1000, ylim = c(-10,10), main = "MDS - Noob - Vape Status -
#  1000 sites")
# 
# mdsPlot(msetNoob, sampGroups = targets$recruitment_center, sampNames = targets$Sample_Name,
#         numPositions=num_1000, ylim = c(-10,10), pal = center_colors, main = "MDS - Noob - Recruitment Center -
#  1000 sites")
# 
# mdsPlot(msetNoob, sampGroups = targets$sex_lab, sampNames = targets$Sample_Name,
#         numPositions=num_10000, ylim = c(-10,10), main = "MDS - Noob - Sex -
#  10000 sites", pal = sex_pal)
# 
# mdsPlot(msetNoob, sampGroups = targets$vape_6mo_lab, sampNames = targets$Sample_Name,
#         numPositions=num_10000, ylim = c(-10,10), main = "MDS - Noob - Vape Status -
#  10000 sites")
# 
# mdsPlot(msetNoob, sampGroups = targets$recruitment_center, sampNames = targets$Sample_Name,
#         numPositions=num_10000, ylim = c(-10,10), pal = center_colors, main = "MDS - Noob - Recruitment Center -
#  10000 sites")
# 
# ################ Normalization and MDS Plots (Autosomes Only) ##################