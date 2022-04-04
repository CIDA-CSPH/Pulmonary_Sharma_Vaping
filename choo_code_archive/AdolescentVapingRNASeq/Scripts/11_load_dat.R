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




metadata_all <-
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



# Possible Vaping Definitions --------------------------------------------------

metadata_all %>%
  select(contains("vape")) %>%
  skim

# vaping in past 30 days (n = 5)
metadata_all$VapingPast30 <-
  ( is.na(metadata_all$Vape_Days) | (metadata_all$Vape_Days == 0) ) %>%
  if_else(., "Control", "Vaping") %>% factor
table(metadata_all$VapingPast30)

# on days vape, how many times used? (n = 18)
metadata_all$VapingTimesPer <-
  is.na(metadata_all$vape_puffs) %>%
  if_else(., "Control", "Vaping") %>% factor
table(metadata_all$VapingTimesPer)

# last vape in last 6 months (n = 12)
metadata_all$VapingLast6Mo <-
  ( is.na(metadata_all$Last_vape) | (metadata_all$Last_vape %in% c(5, 6, 7, 8)) ) %>%
  if_else(., "Control", "Vaping") %>% factor
table(metadata_all$VapingLast6Mo)


#Trent checking Choo's code
metadata_all %>% 
  group_by(VapingLast6Mo) %>% 
  summarise(N = n())








# RNA-seq ======================================================================

# 60651 x 49
count_mat <-
  read_tsv("./Data/20210307_counts/20210307_counts.txt")

# test <-
#   readGFF("Data/20210307_counts/gencode.v37.primary_assembly.annotation.gtf")
# 
# gene_annotations_short <-
#   test %>%
#   transmute(ENSG = gene_id, symbol = gene_name, gene_type) %>%
#   filter(!duplicated(paste0(ENSG, symbol)))
# 
# write_tsv(gene_annotations_short, 
#           "./Data/20210307_counts/gencode_annotations.txt")
gene_annotations_short <- 
  read_tsv("./Data/20210307_counts/gencode_annotations.txt")


# Choo doing some rough gene filtering? (Trent)
# remove counts of 0 (Trent)
filter_genes_percentzero <-
  count_mat[ , -1] %>%
  apply(., 1, function(x) { sum(x == 0) } ) %>%
  `<=`(., 49*0.75)
#remove counts with range > 100? (Trent)
filter_genes_range <-
  count_mat[ , -1] %>%
  apply(., 1, function(x) { max(x) - min(x) } ) %>%
  `>=`(., 100)

filter_genes <-
  filter_genes_percentzero & 
  filter_genes_range



# metadata in order

metadata_naive <-
  read_tsv("./Data/20210307_counts/20201216_coreID_to_PID.txt") %>%
  mutate(NewID = str_pad(NewID, 2, "left", "0") %>% paste0("Sample", .)) %>%
  left_join(metadata_all, by = "SID") %>%
  filter(NewID %in% names(count_mat))

filter_samples <-
  rep(T, length(metadata_naive$SID))



# input into DE ================================================================

#Maybe selecting which samples she is going to use? (Trent)
counts_mat_naive <-
  count_mat[ filter_genes, c(T, filter_samples)]

counts_mat_naive <-
  counts_mat_naive[ , -1] %>%
  as.matrix %>%
  `rownames<-`(., counts_mat_naive$Feature)

metadata_naive <-
  metadata_naive %>%
  filter(filter_samples)

identical(colnames(counts_mat_naive),
          metadata_naive$NewID)

