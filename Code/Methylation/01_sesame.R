####################Load Libraries########################
library(tidyverse)
library(here)
library(kableExtra)
library(sesame) #ver 1.15.6
library(readxl)
library(parallel)
library(randomForest)
library(janitor)

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

#Load EPIC Array (only need to do this the first time)
#sesameDataCache()

#################################### SeSAMe Pipeline By Hand (FINICKY B/C Parallelization!!!!)#######################################
# ################## Read in IDAT data #################
methylation_raw <- lapply(searchIDATprefixes(folder_raw_dat),
                                          readIDATpair)
# 
# ################## Get Metrics  #################
# # methylation_qc <- do.call(rbind, mclapply(methylation_raw, function(x) 
# #                             as_tibble(sesameQC_calcStats(x)))) %>% 
# #   bind_cols(Sample_ID = targets$Sample_Name, .)
# # 
# # rownames(methylation_qc) <- NULL
# 
# ################## Find Outliers  #################
# # outlier_quant <- quantile(methylation_qc$mean_intensity, probs = c(0.01, 0.99))
# # 
# # methylation_qc <- methylation_qc %>% 
# #   mutate(outlier = if_else(mean_intensity < outlier_quant[1], T, F))
# 
# #write_csv(methylation_qc, here("DataProcessed/methylation/methylation_qc_metrics.csv"))
# ################## Plot Outliers  #################
# 
# # methylation_qc %>% 
# #   ggplot(aes(y = log2(mean_intensity), labels = Sample_ID)) +
# #   geom_boxplot()+
# #   ylab("log2(Mean Intensity)")+
# #   theme(axis.title.x=element_blank(),
# #         axis.text.x=element_blank(),
# #         axis.ticks.x=element_blank())
# # 
# # methylation_qc$Sample_ID[methylation_qc$outlier==T]
# # 
# # methylation_qc$mean_intensity[methylation_qc$outlier == T]  
# 
# ################## Dye-Bias Correction Visualization ##############
# par(mfrow=c(2,3))
# sesameQC_plotRedGrnQQ(methylation_raw$`205310770060_R01C01`, main = "R01C01 Before")
# sesameQC_plotRedGrnQQ(methylation_raw$`205310770060_R02C01`, main = "R02C01 Before")
# sesameQC_plotRedGrnQQ(methylation_raw$`205310760169_R03C01`, main = "R03C01 Before")
# sesameQC_plotRedGrnQQ(dyeBiasNL(methylation_raw$`205310770060_R01C01`), main = "R01C01 After")
# sesameQC_plotRedGrnQQ(dyeBiasNL(methylation_raw$`205003700122_R02C01`), main = "R02C01 After")
# sesameQC_plotRedGrnQQ(dyeBiasNL(methylation_raw$`205310760169_R03C01`), main = "R03C01 After")
# 
# ############### Background Subtraction ########################
# par(mfrow=c(2,1), mar=c(3,3,2,1))
# sesameQC_plotBetaByDesign(methylation_raw$`205310770060_R01C01`, main="R01C01 Before", xlab="\beta")
# sesameQC_plotBetaByDesign(noob(methylation_raw$`205310770060_R01C01`), main="R01C01 After", xlab="\beta")
# 
# ################## Run openSesame pipeline step-by-step #################
# #qualityMask
# methylation_mask <- BiocParallel::bplapply(methylation_raw, qualityMask)
# #Infer Infinium I Channel
# methylation_infer <- BiocParallel::bplapply(methylation_raw, inferInfiniumIChannel)
# #Dye-Bias Correction
# methylation_db_corr <- BiocParallel::bplapply(methylation_infer, dyeBiasNL)
# #detection p-value masking
# methylation_detect_mask <- BiocParallel::bplapply(methylation_db_corr, pOOBAH)
# #Background subtraction
# methylation_background <- BiocParallel::bplapply(methylation_detect_mask, noob)
# #getBetas
# betas <- do.call(cbind, BiocParallel::bplapply(methylation_background, getBetas))
# 
# betas <- openSesame(searchIDATprefixes(folder_raw_dat), BPPARAM = BiocParallel::MulticoreParam(2))
# 
# betas <- betas %>% 
#   data.frame %>% 
#   mutate(CpG_Site = rownames(.)) %>% 
#   select(CpG_Site, everything())

################# Check Probe Matching #####################
par(mfrow=c(2,1), mar=c(3,3,2,1))
sesameQC_plotBetaByDesign(methylation_raw$`205310770060_R01C01`, main="Before", xlab="\beta")
sesameQC_plotBetaByDesign(matchDesign(methylation_raw$`205310770060_R01C01`), main="After", xlab="\beta")

################# Write out Betas #####################

#write_tsv(betas, file = here("DataProcessed/methylation/methylation_betas.txt"))

#fwrite(betas, file = here("DataProcessed/methylation/methylation_betas_full.txt"))


############## Convert to M Values #####################
mvals <- betas %>% 
  select(-CpG_Site) %>% 
  apply(., 2, BetaValueToMValue) %>% 
  as.data.frame() %>% 
  mutate(CpG_Site = rownames(.)) %>% 
  select(CpG_Site, everything())

#fwrite(mvals, here("DataProcessed/methylation/methylation_mvals.txt"))