# working directory & libraries ------------------------------------------------

date_export <-
  "20210326"

setwd("/beevol/home/borengas/analyses/methylation")

library(tidyverse)
library(readxl)
library(ggthemes)
library(ggbeeswarm)
library(cowplot)

library(data.table)
library(sesameData)
library(cowplot)
library(ggrepel)
library(ggthemes)

library(qvalue)
library(vegan)
library(RUVSeq)

library(bacon)
library(lmerTest)



# palettes  --------------------------------------------------------------------




# palettes =====================================================================

palette_color_vaping <-
  c(Control = "black", Vaping = "#b2182b")
palette_color_center <-
  c(Pueblo = "#e41a1c", `CommCity/Denver` = "#377eb8", Aurora = "#4daf4a")
palette_shape_vaping <-
  c(Control = 1, Vaping = 19)









# metadata =====================================================================
# from Daniel
# one sample ID = "DOO3" non-numeric
dat_sample_qc <- 
  read_excel("./Data/20201208/RNA-DNA extraction Vaping 11_23_2020.xlsx")

# corrected document, with concentrations
dat_qc_RNA <-
  read_excel("./Data/20201217/RNA numbers_Vaping Samples SS.xlsx") %>%
  dplyr::slice(-1) %>%
  select(-1) %>%
  transmute(SID = `...3`,
            `260/280` = `RNA...4`,
            `Conc (ng/uL)` = `RNA...5`,
            `Total RNA (ng)` = `Total RNA`) %>%
  mutate_all(as.double) %>% # one NA, "DOO3"
  filter(!is.na(`260/280`) & !is.na(SID))

dat_participants <-
  read_csv("./Data/20210325/all_data_Latino_Youth_survey_clean_pulmfxn_FINAL_3.23.21.csv")



final_metadata <-
  dat_participants %>%
  left_join(., dat_qc_RNA) %>%
  mutate(CityAurora = grepl("Aurora", City),
         CityPueblo = grepl("Pueblo", City),
         RecruitmentCenter = case_when(
           grepl("Aurora", City) ~ "Aurora",
           grepl("Pueblo|Avondale|Mineral", City) ~ "Pueblo",
           T ~ "CommCity/Denver"
         ),
         GenderMale = factor(Gender == "1"),
         EthnicityLatino = factor(Latino == "1"),
         X20X5 = X20 - X5,
         R20R5 = R20 - R5)



final_metadata <-
  read_tsv("./Data/20210401_methylation_coreID_to_PID.txt") %>%
  mutate(NewID = str_pad(NewID, 2, "left", "0") %>% paste0("Sample", .)) %>%
  left_join(final_metadata, by = "SID")





# Possible Vaping Definitions --------------------------------------------------

final_metadata %>%
  select(contains("vape")) %>%
  skim

# vaping in past 30 days (n = 5)
final_metadata$VapingPast30 <-
  ( is.na(final_metadata$Vape_Days) | (final_metadata$Vape_Days == 0) ) %>%
  if_else(., "Control", "Vaping") %>% factor
table(final_metadata$VapingPast30)

# on days vape, how many times used? (n = 18)
final_metadata$VapingTimesPer <-
  is.na(final_metadata$vape_puffs) %>%
  if_else(., "Control", "Vaping") %>% factor
table(final_metadata$VapingTimesPer)

# last vape in last 6 months (n = 12)
final_metadata$VapingLast6Mo <-
  ( is.na(final_metadata$Last_vape) | (final_metadata$Last_vape %in% c(5, 6, 7, 8)) ) %>%
  if_else(., "Control", "Vaping") %>% factor
table(final_metadata$VapingLast6Mo)


