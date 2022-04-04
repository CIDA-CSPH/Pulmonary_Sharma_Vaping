# Possible Vaping Definitions --------------------------------------------------

metadata <-
  dat_participants %>%
  left_join(., dat_qc_DNA) %>%
  filter(!is.na(`Total DNA (ng)`)) %>%
  filter(SID != "102") # tech issues/putative outlier in RNA-seq data

# last vape in last 6 months (n = 12)
metadata$VapingLast6Mo <-
  ( is.na(metadata$Last_vape) | (metadata$Last_vape %in% c(5, 6, 7, 8)) ) %>%
  if_else(., "Control", "Vaping")
table(metadata$VapingLast6Mo)




rand_Age <-
  ( metadata$Age <= median(metadata$Age) ) %>%
  if_else(., "AgeLow", "AgeHigh")
rand_gender <-
  (metadata$Gender == 1) %>%
  if_else(., "M", "F")
rand_ethn <-
  metadata$Latino
rand_vape <-
  metadata$VapingLast6Mo

randgroup <-
  # paste0(rand_Age, rand_vape) %>%
  paste0(rand_Age, rand_ethn, rand_gender, rand_vape) %>%
  as.factor
table(randgroup)

randgroup_num <-
  randgroup %>%
  set_fctlevels(1:14)


ggplot(data = NULL,
       aes(x = 1:50, y = as.numeric(randgroup_num))) +
  geom_point() +
  geom_line(group = 1)

ggplot(data = NULL,
       aes(x = 1:50, y = metadata$Age + 7)) +
  geom_point() +
  geom_line(group = 1)


library(ggridges)
ggplot(data = metadata,
       aes(Age + 7)) +
  geom_bar() +
  facet_grid(. ~ (City == "Pueblo")) +
  scale_x_continuous(breaks = 5:10 + 7)
  # geom_density_ridges(
  #   jittered_points = TRUE,
  #   position = position_points_jitter(width = 0.05, height = 0),
  #   point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
  # )

library(ggridges)
ggplot(data = metadata,
       aes(SID, City == "Pueblo")) +
geom_density_ridges(
  jittered_points = TRUE,
  position = position_points_jitter(width = 0.05, height = 0),
  point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
)

table(metadata$Sex, metadata$VapingLast6Mo)

levels(randgroup)[c(5, 7)]










# randomize by 6 (hypothetical for methylation 8 slots per chip, ultimately)

tmp <- 1:nrow(metadata) # should be 1:48
randomize_by_six_check_seed <-
  function(seed, check = T) {
    
    set.seed(seed)
    randomization_chip <-
      createFolds(randgroup, 6)
  
    sapply(
      1:length(randomization_chip),
      function(x) {
      tmp <- get("tmp", envir = .GlobalEnv)
      tmp[ randomization_chip[[x]] ] <- paste0("Chip", x)
      assign("tmp", tmp, envir = .GlobalEnv)
      ""
  })
  
    N_vapers_by_batch <-
      metadata %>%
      bind_cols(Chip = tmp) %>%
      filter(VapingLast6Mo == "Vaping") %>%
      group_by(Chip) %>%
      summarize(N = n()) %>%
      .$N
    
    
    N_by_batch <-
      metadata %>%
      bind_cols(Chip = tmp) %>%
      group_by(Chip) %>%
      summarize(N = n()) %>%
      .$N
  
    all(N_vapers_by_batch %in% 1:3) &
      all(N_by_batch == 8)  
    
}

# (
# valid_seeds <-
#   sapply(1:1000, randomize_by_six_check_seed) %>%
#   which )

# valid_seeds
randomize_by_six_check_seed(194)
tmp %>% table
final <-
  metadata %>%
  bind_cols(Chip = tmp) 





randomize_by_column_check_seed <- function(seed) {
  
  set.seed(seed)
  
  final_ordering <-
    final %>%
    group_by(Chip) %>%
    group_split() %>%
    map(~ bind_cols(.x, RowWithinChip = rnorm(n = nrow(.x)))) %>%
    map(~ mutate(.x, RowWithinChip = order(RowWithinChip))) %>%
    bind_rows() %>%
    arrange(Chip, RowWithinChip)
  
  assign("final_ordering", final_ordering, envir = .GlobalEnv)
  
  final_ordering %>%
    group_by(RowWithinChip) %>%
    summarize(N = n(),
              sum(City == "Pueblo")/N,
              sum(Latino == "1")/N,
              sum(Gender == "1")/N,
              sum(VapingLast6Mo == "Vaping")/N) %>%
    select(-(1:2)) %>%
    mutate_all(function(x) { x != 0 & x != 1 } ) %>%
    unlist %>%
    all
}

# ( seeds <-
#   sapply(1:500, randomize_by_column_check_seed) %>%
#   which )

randomize_by_column_check_seed(4)






final_ordering %>%
  group_by(Chip) %>%
  summarize(N = n(),
            `Age (yrs)` = summaryMedianIQR(Age),
            `City (Pueblo)` = summaryCountPercent(City, "Pueblo"),
            `Ethnicity (Latino)` = summaryCountPercent(Latino, "1"),
            `Gender (Male)` = summaryCountPercent(`Gender`, "1"),
            `Vaping (last 6mo)` = summaryCountPercent(VapingLast6Mo, "Vaping")) %>%
  View

final_ordering %>%
  group_by(RowWithinChip) %>%
  summarize(N = n(),
            `Age (yrs)` = summaryMedianIQR(Age),
            `City (Pueblo)` = summaryCountPercent(City, "Pueblo"),
            `Ethnicity (Latino)` = summaryCountPercent(Latino, "1"),
            `Gender (Male)` = summaryCountPercent(`Gender`, "1"),
            `Vaping (last 6mo)` = summaryCountPercent(VapingLast6Mo, "Vaping")) %>%
  View()

for_daniel <-
  final_ordering %>%
  arrange(Chip, RowWithinChip) %>%
  transmute(SID,
            NewID = 1:nrow(.),
         `96Well_Row` = rep(LETTERS[1:8], times = 7)[1:48],
         `96Well_Column` = gsub("Chip", "", Chip))

for_daniel_wide <-
  for_daniel %>%
  mutate(info = paste0("SID = ", SID, "\nNewID = ", NewID)) %>%
  pivot_wider(id_cols = `96Well_Row`,
              names_from = `96Well_Column`,
              names_prefix = "Column_",
              values_from = info)

for_daniel %>%
  transmute(NewID, SID) %>%
  write_tsv("./Output/20210401_methylation_coreID_to_PID.txt")


for_core <-
  for_daniel %>%
  left_join(dat_qc_DNA %>% select(SID, `Total DNA (ng)`)) %>%
  select(NewID, `96Well_Row`, `96Well_Column`, `Total DNA (ng)`) %>%
  set_names(c("ID", "Row", "Column/Chip", "Total DNA (ng)"))


for_core_wide <-
  for_core %>%
  pivot_wider(id_cols = `Row`,
              names_from = `Column/Chip`, names_prefix = "Chip_",
              values_from = ID)
  
View(for_core)
View(for_core_wide)

for_core$ID
  






write_xlsx(list(SharmaVapingMethylation = for_daniel,
                wide = for_daniel_wide),
           path = "Output/20210401-SharmaVapingMethylation-forDW.xlsx")

write_xlsx(list(SharmaVapingMethylation = for_core,
                wide = for_core_wide),
           path = "Output/20210401-SharmaVapingMethylation-forcore.xlsx")


