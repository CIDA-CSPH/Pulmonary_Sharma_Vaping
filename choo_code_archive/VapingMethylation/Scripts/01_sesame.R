# ==============================================================================
# 01_sesame_pipeline.R
# .idat --> beta values
# ==============================================================================



# server setup -----------------------------------------------------------------

screen
qlogin -R rusage[mem=48]
cd /beevol/home/lincuini/analyses/
  module load gcc/7.4.0
module load R/4.0.3
R



# libraries --------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(readxl)
library(data.table)
library(Cairo)
library(minfi)
library(sesame)
library(khroma)


# sample metadata --------------------------------------------------------------
setwd("/beevol/home/lincuini/Analyses/VapingMethylation/")
date_export <- "20210914"
folder_raw_dat <-
  "./Data/Idats/"

samplesheet <-
  left_join(
    read_csv("./Data/Metadata/SSharma_48_samplesheet_05202021.csv", skip = 7) %>%
      select(-`Sample_Group`, -`Pool_ID`),
    read_excel("./Data/Metadata/Sharma_48_5192021Controls.xlsx") %>%
      mutate(`Sample Name` = as.numeric(`Sample Name`)),
    by = c("Sample_Name" = "Sample Name"))
# 
# 
# # manual entry from Okyong Cho e-mail
# # Thursday, October 8, 2020 3:39 PM
# assay_scan_date <-
#   tibble(EPIC_Date_Scanned = c(rep("09042019", 6), rep("09052019", 4), rep("09062019", 6)),
#          EPIC_chip = c("203784950100", "203784950110", "203806760084", "203806760086", "203806760118", "203806760124",
#                        "203784950152", "203806760090", "203806760119", "203806760153",
#                        "203784950101", "203784950111", "203784950132", "203784950146", "203784950151", "203806760151"))


# file path example formats
# Slide_1 _203784950100/203784950100_R01C01_Grn.idat

targets <-
  tibble(file_path =
           list.files(path = folder_raw_dat, pattern = "*.idat|*.IDAT",
                      full.names = F, recursive = T)) %>%
  mutate(`Sentrix Barcode` = str_split(file_path, pattern = "/|_") %>% map_chr( ~ .[4]), # chip
         `Sentrix Position` = str_split(file_path, pattern = "/|_") %>% map_chr( ~ .[5])) %>% # row
  left_join(samplesheet,
            by = c("Sentrix Barcode", "Sentrix Position")) %>%
  mutate(file_path_prefix = gsub("_Red.idat|_Grn.idat", "", file_path) %>%
           paste0(folder_raw_dat, .)) %>%
  filter(!duplicated(file_path_prefix)) #%>% # duplicated prefix due to Red/Grn
  # left_join(assay_scan_date, by = "EPIC_chip") %>%
  # mutate(EPIC_WellPairs = EPIC_row %>% as.factor() %>%
  #          fct_collapse(`1/2` = c("R01C01", "R02C01"),
  #                       `3/4` = c("R03C01", "R04C01"),
  #                       `5/6` = c("R05C01", "R06C01"),
  #                       `7/8` = c("R07C01", "R08C01")))
  # 



# sesameDataCache("EPIC")



methylation_raw <- 
  map(targets$file_path_prefix, readIDATpair)

sesame_qc <-
  methylation_raw %>%
  map(sesameQC) %>%
  map_dfr( ~ as_tibble(.x)) %>%
  bind_cols(ID = targets$Sample_Name, .)

methylation_betas <-
  targets$file_path_prefix %>%
  openSesame()

probenames <- rownames(methylation_betas)

methylation_betas <-
  methylation_betas %>% 
  as.data.table %>%
  set_names(., targets$Sample_Name %>% paste0("Sample", .))

eps <- 1e-6
methylation_mvals <- copy(methylation_betas)
methylation_mvals <-
  methylation_mvals %>%
  .[, (names(.)) := lapply(.SD, function(x) { log( (x + eps) / ( 1 - x + eps) ) } )]


dir.create("Output")
dir.create(paste0("Output/", date_export, "_MethylationDat"))

fwrite(
  cbind(Probe = probenames, methylation_betas),
  file = paste0("./Output/", date_export, "_MethylationDat/betas.tsv"),
  sep = "\t")
fwrite(
  cbind(Probe = probenames, methylation_mvals),
  file = paste0("./Output/", date_export, "_MethylationDat/mvals.tsv"),
  sep = "\t")

write_tsv(
  sesame_qc,
  path = paste0("./Output/", date_export, "_MethylationDat/sesame_qc.tsv"))
write_tsv(
  targets,
  path = paste0("./Output/", date_export, "_MethylationDat/targets.tsv"))
methylation_raw %>% map_dfc(~ .x@ctl[ , "G"]) %>% write_tsv(
  path = paste0("./Output/", date_export, "_MethylationDat/controls_G.tsv")
)
methylation_raw %>% map_dfc(~ .x@ctl[ , "R"]) %>% write_tsv(
  path = paste0("./Output/", date_export, "_MethylationDat/controls_R.tsv")
)

save(
  methylation_raw, methylation_betas, targets, sesame_qc,
  file = paste0("./Output/", date_export, "_MethylationDat/sesameoutput.Rdata"))







# 8 samples per plate, 16 plates  
palette_colors <-
  c(colour("light")(6))#, colour("muted")(8))

export_expression_boxplots <-
  function(dat, fileout, ylabel) {
    par(mar = c(.2, .2, .2, .2))
    CairoPNG(filename = paste0("./Output/", date_export, "_MethylationDat/", fileout),
             dpi = 300, width = 1600, height = 1600)
    range <- dat %>% range(na.rm = T)
    layout(matrix(c(1, 1, 1,
                    2, 2, 2),
                  nrow = 2, byrow = TRUE))
    boxplot(dat[ , 1:24],
            boxfill = rep(palette_colors[1:3], each = 8),
            xaxt = "n", ylim = range, ylab = ylabel)
    boxplot(dat[ , 25:48],
            boxfill = rep(palette_colors[4:6], each = 8),
            xaxt = "n", ylim = range, ylab = ylabel)
    dev.off()
  }

export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@II[ , "M"]),
                           "rawintensities_II_M.png", "Type II (M)")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@II[ , "U"]),
                           "rawintensities_II_U.png", "Type II (U)")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@IG[ , "M"]),
                           "rawintensities_IG_M.png", "Type IG (M)")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@IG[ , "U"]),
                           "rawintensities_IG_U.png", "Type IG (U)")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@oobG[ , "M"]),
                           "rawintensities_oobG_M.png", "Type OOBG (M)")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@oobG[ , "U"]),
                           "rawintensities_oobG_U.png", "Type OOBG (U)")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@oobR[ , "M"]),
                           "rawintensities_oobR_M.png", "Type OOBR (M)")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@oobR[ , "U"]),
                           "rawintensities_oobR_U.png", "Type OOB (U)")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@ctl[ , "G"]),
                           "rawintensities_ctl_G.png", "Control Green")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x@ctl[ , "R"]),
                           "rawintensities_ctl_R.png", "Control Red")


set.seed(1234)
subsample_index <- sample(nrow(methylation_betas), 10000)
export_expression_boxplots(
  methylation_betas[subsample_index, ] %>% as.matrix,
  "processed_betas.png",
  "SeSAMe-Normalized Betas")
export_expression_boxplots(
  methylation_mvals[subsample_index, ] %>% as.matrix,
  "processed_mvals.png",
  "SeSAMe-Normalized M-Values")

CairoPNG(filename = paste0("./Output/", date_export, "_MethylationDat/example_densities_betas.png"),
         dpi = 300, width = 1000, height = 1000)
par(cex = 0.5)
densityPlot(methylation_betas %>% as.matrix,
            sampGroups = targets$`Sentrix Barcode`,
            pal = palette_colors)
dev.off()

CairoPNG(filename = paste0("./Output/", date_export, "_MethylationDat/example_densities_mvals.png"),
         dpi = 300, width = 1000, height = 1000)
par(cex = 0.5)
densityPlot(methylation_mvals %>% as.matrix,
            sampGroups = targets$`Sentrix Barcode`,
            pal = palette_colors)
dev.off()





# quit()
# n
# exit


