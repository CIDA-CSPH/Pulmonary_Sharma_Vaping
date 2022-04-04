# from Daniel
# one sample ID = "DOO3" non-numeric
dat_sample_qc <- 
  read_excel("./Data/20201208/RNA-DNA extraction Vaping 11_23_2020.xlsx")

# isolate DNA-only component
dat_qc_DNA <-
  dat_sample_qc %>%
  slice(-1) %>%
  transmute(SID = `...8`,
            `260/280` = DNA,
            `Total DNA (ng)` = `Total DNA`) %>%
  mutate_all(as.double) %>% # one NA, "DOO3"
  filter(!is.na(`260/280`) & !is.na(SID)) # 0 ng, empty vial

# RNA-only
# corrected document, has fixed concentrations
dat_qc_RNA <-
  read_excel("./Data/20201217/RNA numbers_Vaping Samples SS.xlsx") %>%
  slice(-1) %>%
  select(-1) %>%
  transmute(SID = `...3`,
            `260/280` = `RNA...4`,
            `Conc (ng/uL)` = `RNA...5`,
            `Total RNA (ng)` = `Total RNA`) %>%
  mutate_all(as.double) %>% # one NA, "DOO3" from different project
  filter(!is.na(`260/280`) & !is.na(SID))


dat_participants <-
  read_csv("./Data/20201208/all_data_Latino_Youth_survey_clean_pulmfxn_FINAL.csv")


metadata <-
  dat_participants %>%
  left_join(., dat_qc_RNA) 



