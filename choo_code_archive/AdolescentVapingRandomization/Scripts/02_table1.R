# Possible Vaping Definitions --------------------------------------------------

metadata %>%
  select(contains("vape")) %>%
  skim

# vaping in past 30 days (n = 5)
metadata$VapingPast30 <-
  ( is.na(metadata$Vape_Days) | (metadata$Vape_Days == 0) ) %>%
  if_else(., "Control", "Vaping")
table(metadata$VapingPast30)

# on days vape, how many times used? (n = 18)
metadata$VapingTimesPer <-
  is.na(metadata$vape_puffs) %>%
  if_else(., "Control", "Vaping")
table(metadata$VapingTimesPer)

# last vape in last 6 months (n = 12)
metadata$VapingLast6Mo <-
  ( is.na(metadata$Last_vape) | (metadata$Last_vape %in% c(5, 6, 7, 8)) ) %>%
  if_else(., "Control", "Vaping")
table(metadata$VapingLast6Mo)



# table 1 ----------------------------------------------------------------------



# metadata$vapinggroup <-
#   metadata$VapingTimesPer



# IOS = https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4948229/
# age, sex dependent?
# - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3401927/
# - https://www.intechopen.com/books/asthma-biological-evidences/clinical-applications-of-impulse-oscillometry

metadata$vapinggroup <-
  metadata$VapingLast6Mo

metadata %>%
  group_by(vapinggroup) %>%
  mutate(`R20 - R5` = R20 - R5,
         `X20 - X5` = X20 - X5) %>%
  summarise(
    
    `Basic Demographics` = summaryBlankLine(),
    
    N = n() %>% as.character(),
    
    `Age (years)` = summaryMeanSD(Age + 7),
    `Gender (Male)` = summaryCountPercent(Gender, "1"),
    `Ethnicity (Hispanic/Latino)` = summaryCountPercent(Latino, "1"),
    # `Race (American-Indian/Alaska Native)` = summaryCountPercent(Race, "1"),
    `Race (Asian)` = summaryCountPercent(Race, "2"),
    `  (African-American/Black)` = summaryCountPercent(Race, "3"),
    `  (Native Hawaiian/Pacific Islander)` = summaryCountPercent(Race, "4"),
    `  (White)` = summaryCountPercent(Race, "5"),
    `  (Other)` = summaryCountPercent(Race, "7"),
    
    `City (Aurora)` = summaryCountPercent(City, "Aurora", fuzzy = T),
    `   (Commerce City)` = summaryCountPercent(City, "Commerce", fuzzy = T),
    `   (Denver)` = summaryCountPercent(City, "Denver", fuzzy = T),
    `   (Pueblo)` = summaryCountPercent(City, "Pueblo", fuzzy = T),
    
  
    `Vaping/Smoke Exposure` = summaryBlankLine(),
    
    `Age Began Vaping (years)` = summaryMeanSD(Vape_firstage, na.rm = T),
    `Ever Vaped` = summaryCountPercent(Ever_vape, 1),
    `One Days Vape, Times Pickup Vape / Day` = summaryMeanSD(Vape_pickup, na.rm = T),
    `Each Time Vape, Number Puffs / Use` = summaryMeanSD(Vape_puffs, na.rm = T),
    `Frequency Refill (At Least Once/Day)` = summaryCountPercent(Freq_refill, "1"),
    `   (Every Couple Days or Once A Week)` = summaryCountPercent(Freq_refill, c("2", "3")),
    `Any Cigarettes / Past 30 Days` = summaryCountPercent(cig_30days, c("7", "8"), inverse = T),
    `Any Use / Past 30 Days` = summaryCountPercent(mj_30days, c("2", "8"), inverse = T),
    `Vape Contains Nicotine` = summaryCountPercent(Nicotine, "1"),
    `Ever Tobacco` = summaryCountPercent(ever_tobacco, "10", inverse = T),
    `Ever Cigarettes` = summaryCountPercent(cig_life, "1"),
    `Ever Vape` = summaryCountPercent(Ever_vape, "1"),
    
    
    `Spirometry/Oscillometry Measurements` = summaryBlankLine(),
    `FEV1(pp%)` = summaryMeanSD(`FEV1%`, na.rm = T),
    `FVC(pp%)` = summaryMeanSD(`FVC%`, na.rm = T),
    `FEV1/FVC` = summaryMeanSD(`FEV1/FVC`, na.rm = T, digits = 2),
    `Exhaled NO` = summaryMeanSD(`Exhaled NO`, na.rm = T, digits = 2),
    `R5` = summaryMeanSD(R5, na.rm = T, digits = 2),
    `R20` = summaryMeanSD(R20, na.rm = T, digits = 2),
    `R20 - R5` = summaryMeanSD(`R20 - R5`, na.rm = T, digits = 2),
    `X5` = summaryMeanSD(X5, na.rm = T, digits = 2),
    `X20` = summaryMeanSD(X20, na.rm = T, digits = 2),
    `X20 - X5` = summaryMeanSD(`X20 - X5`, na.rm = T, digits = 2),
    `Tidal Volume (VT)` = summaryMeanSD(VT, na.rm = T, digits = 2),
    `Area of Reactance (AX)` = summaryMeanSD(AX, na.rm = T, digits = 2),
    `Resonant Frequency (fres)` = summaryMeanSD(Fres, na.rm = T, digits = 2),
    
    `Respiratory Symptoms` = summaryBlankLine(),
    `History Lung Disease` = summaryCountPercent(lung_dis, "1"),
    `Cough in Last 2 Weeks` = summaryCountPercent(cough, "1"),
    `Ever Wheezing While Breathing` = summaryCountPercent(wheezing, "1"),
    `Shortness Breath in Last 2 Weeks` = summaryCountPercent(SOB, "1"),
    `Ever Difficult Breathing, Pressure On Chest` = summaryCountPercent(breathe_diff, "1"),
  ) %>%
  pivot_longer(cols = 2:ncol(.), names_to = " ") %>%
  pivot_wider(id_cols = ` `, names_from = vapinggroup)  %>%
  mutate(` ` = if_else(Vaping != "", paste0(". . . ", ` `), ` `)) %>%
  
  View(title = "Table 1")

