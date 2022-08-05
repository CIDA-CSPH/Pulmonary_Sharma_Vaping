####################Load Libraries########################
library(tidyverse)
library(here)
library(kableExtra)
library(sesame)
library(readxl)
library(parallel)

kablize <- function(tab, digits = 3) {
  kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}

format_num <- function(number, digits = 0, format = "f") {
  formatC(number, digits = digits, format = format, big.mark = ",")
}
################## Setup #################
#Read in and join metadata
samplesheet <- left_join(
    read_csv(here("DataRaw/methylation/Metadata/SSharma_48_samplesheet_05202021.csv"), skip = 7) %>%
      select(-`Sample_Group`, -`Pool_ID`),
    read_excel(here("DataRaw/methylation/Metadata/Sharma_48_5192021Controls.xlsx")) %>%
      mutate(`Sample Name` = as.numeric(`Sample Name`)),
    by = c("Sample_Name" = "Sample Name"))

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

#Load EPIC Array
sesameDataCache("EPIC")

################## Read in IDAT data #################
methylation_raw <- mclapply(searchIDATprefixes(folder_raw_dat), 
                                          readIDATpair, 
                                          mc.cores = 2)

################## Get Metrics  #################
methylation_qc <- do.call(rbind, lapply(methylation_raw, function(x)
  as_tibble(sesameQC_calcStats(x)))) %>% 
  bind_cols(Sample_ID = targets$Sample_Name, .)

rownames(methylation_qc) <- NULL

################## Find Outliers  #################
outlier_quant <- quantile(methylation_qc$mean_intensity, probs = c(0.01, 0.99))

methylation_qc <- methylation_qc %>% 
  mutate(outlier = if_else(mean_intensity < outlier_quant[1], T, F))

################## Plot Outliers  #################

methylation_qc %>% 
  ggplot(aes(y = log2(mean_intensity), labels = Sample_ID)) +
  geom_boxplot()+
  ylab("log2(Mean Intensity)")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

methylation_qc$Sample_ID[methylation_qc$outlier==T]

methylation_qc$mean_intensity[methylation_qc$outlier == T]  

################## Run openSesame pipeline step-by-step #################
#qualityMask
methylation_mask <- BiocParallel::bplapply(methylation_raw, qualityMask)
#Infer Infinium I Channel
methylation_infer <- BiocParallel::bplapply(methylation_mask, inferInfiniumIChannel)
#Dye-Bias Correction
methylation_db_corr <- BiocParallel::bplapply(methylation_infer, dyeBiasNL)
#detection p-value masking
methylation_detect_mask <- BiocParallel::bplapply(methylation_db_corr, pOOBAH)
#Background subtraction
methylation_background <- BiocParallel::bplapply(methylation_detect_mask, noob)
#getBetas
betas <- do.call(cbind, BiocParallel::bplapply(methylation_background, getBetas))

################## Compare to OpenSesame #################
betas_test <- openSesame(folder_raw_dat, BPPARAM = BiocParallel::MulticoreParam(2))
#Already ran the above line. Results match exactly. 

################## Dye-Bias Correction Visualization ##############

sesamePlotRedGrnQQ(methylation_raw$`205310770060_R01C01`)
sesamePlotRedGrnQQ(methylation_db_corr$`205310770060_R01C01`)

par(mfrow=c(1,2))
sesamePlotRedGrnQQ(dyeBiasCorr(methylation_raw$`205310770060_R01C01`)) # linear correction
sesamePlotRedGrnQQ(dyeBiasNL(methylation_raw$`205310770060_R01C01`))   # nonlinear correction

############### Infinium Matching #######################
