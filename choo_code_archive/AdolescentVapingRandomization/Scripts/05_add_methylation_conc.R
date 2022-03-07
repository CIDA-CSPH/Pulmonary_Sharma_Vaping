library(tidyverse)
library(readxl)
library(writexl)
concs <-
  read_excel("./Data/20210416/DNA Extraction_Vaping Project_04-16-21.xlsx",
             skip = 2)



forDW_addconc <-
  read_excel("./Output/20210401-SharmaVapingMethylation-forDW.xlsx",
             skip = 10) %>%
  left_join(., concs %>% mutate(ID = as.double(ID)), by = c("SID" = "ID"))


forcore_addconc <-
  read_excel("./Output/20210401-SharmaVapingMethylation-forcore.xlsx") %>%
  left_join(., forDW_addconc %>% transmute(NewID, `Conc (ng/uL)` = `ng/uL`),
            by = c("ID" = "NewID")) %>%
  select(-`Total DNA (ng)`)




write_xlsx(list(Samples = forcore_addconc,
                wide = read_excel("./Output/20210401-SharmaVapingMethylation-forcore.xlsx", sheet = 2)),
            "./Output/20210416-SharmaVapingMethylation-forcore.xlsx")
