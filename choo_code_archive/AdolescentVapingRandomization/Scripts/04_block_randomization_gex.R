metadata <-
  metadata %>%
  filter(!is.na(`Total RNA (ng)`))

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

Nsamples <- length(randgroup)
ggplot(data = NULL,
       aes(x = 1:Nsamples, y = as.numeric(randgroup_num))) +
  geom_point() +
  geom_line(group = 1)

ggplot(data = NULL,
       aes(x = 1:Nsamples, y = metadata$Age + 7)) +
  geom_point() +
  geom_line(group = 1)

# note that pueblo tends to be older, tend to be vapers
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

tmp <- double(length = Nsamples)
randomize_by_six_check_seed <-
  function(seed, check = T) {
    
    
    set.seed(seed)
    randomization_batch <-
      createFolds(randgroup, 6)
  
    sapply(
      1:length(randomization_batch),
      function(x) {
      tmp <- get("tmp", envir = .GlobalEnv)
      tmp[ randomization_batch[[x]] ] <- paste0("Batch", x)
      assign("tmp", tmp, envir = .GlobalEnv)
      ""
  })
  
    N_vapers_by_batch <-
      metadata %>%
      bind_cols(Batch = tmp) %>%
      filter(VapingLast6Mo == "Vaping") %>%
      group_by(Batch) %>%
      summarize(N = n()) %>%
      .$N
    
    
    N_by_batch <-
      metadata %>%
      bind_cols(Batch = tmp) %>%
      group_by(Batch) %>%
      summarize(N = n()) %>%
      .$N
  
    all(N_vapers_by_batch %in% 1:3) &
      all(N_by_batch %in% 8:9)
}

# (
# valid_seeds <-
#   sapply(1:1000, randomize_by_six_check_seed) %>%
#   which )

# valid_seeds - choose 576

# 194 576 735 963 977 979


randomize_by_six_check_seed(576)
final <-
  metadata %>%
  bind_cols(Batch = tmp) 

check_distribution_by_batch()



  
check_distribution_by_batch <- function() {
final %>%
  group_by(Batch) %>%
  summarize(N = n(),
            `Age (yrs)` = summaryMedianIQR(Age),
            `City (Pueblo)` = summaryCountPercent(City, "Pueblo"),
            `Ethnicity (Latino)` = summaryCountPercent(Latino, "1"),
            `Gender (Male)` = summaryCountPercent(`Gender`, "1"),
            `Vaping (last 6mo)` = summaryCountPercent(VapingLast6Mo, "Vaping")) %>%
  View
}



check_distribution_by_batch()

# Batch4 has 9 samples instead of 8
# move 145, non-Latino, Commerce City
final %>%
  filter(Batch %in% c("Batch4")) %>%
  arrange(Batch, VapingLast6Mo, Age) %>%
  select(Batch, SID, VapingLast6Mo, Age, City, Latino, Gender)

final <-
  final %>% 
  mutate(Batch = if_else(SID %in% c(145), "Batch7", Batch))
  
final %>%
  filter(Batch == "Batch7") %>%
  select(Batch, SID, VapingLast6Mo, Age, City, Latino, Gender)
# if plated by row instead of column

  
  
check_distribution_by_batch()
  





# now check order within batch so vapers aren't next to each other
# randomize withinbatch
# seed = 3 seems to reasonably randomly distribute all variables across col
# manual criteria: at least one vape-exposed in each column, no 100% of any category
randomize_by_column_check_seed <- function(seed) {

  set.seed(seed)
  
  final_ordering <-
    final %>%
    group_by(Batch) %>%
    group_split() %>%
    map(~ bind_cols(.x, ColumnWithinBatch = rnorm(n = nrow(.x)))) %>%
    map(~ mutate(.x, ColumnWithinBatch = order(ColumnWithinBatch))) %>%
    bind_rows() %>%
    arrange(Batch, ColumnWithinBatch)

  assign("final_ordering", final_ordering, envir = .GlobalEnv)
  
  final_ordering %>%
    group_by(ColumnWithinBatch) %>%
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

# seeds <-
#   sapply(1:500, randomize_by_column_check_seed) %>%
#   which
# seeds
# [1]  26  77  92 139 156 169 222 230 242 246 304 332 365 367 397 420 486 495

randomize_by_column_check_seed(26)

final_ordering %>%
  group_by(ColumnWithinBatch) %>%
  summarize(N = n(),
            `Age (yrs)` = summaryMedianIQR(Age),
            `City (Pueblo)` = summaryCountPercent(City, "Pueblo"),
            `Ethnicity (Latino)` = summaryCountPercent(Latino, "1"),
            `Gender (Male)` = summaryCountPercent(`Gender`, "1"),
            `Vaping (last 6mo)` = summaryCountPercent(VapingLast6Mo, "Vaping")) %>%
  View

final_ordering %>%
  group_by(Batch) %>%
  summarize(N = n(),
            `Age (yrs)` = summaryMedianIQR(Age),
            `City (Pueblo)` = summaryCountPercent(City, "Pueblo"),
            `Ethnicity (Latino)` = summaryCountPercent(Latino, "1"),
            `Gender (Male)` = summaryCountPercent(`Gender`, "1"),
            `Vaping (last 6mo)` = summaryCountPercent(VapingLast6Mo, "Vaping")) %>%
  View


old_incorrect_IL5 <-
  read_excel("Output/20201215_Vaping_RNA_Randomize.xlsx") %>%
  mutate(`Bad Relabel (Mistake)` = if_else(NewID %in% 1:20, paste0(NewID, "*"), "") %>%
           if_else(`%in%`(., paste0(c(1:9), "*")), paste0(0, .), .))


for_daniel <-
  final_ordering %>%
  arrange(Batch, ColumnWithinBatch) %>%
  mutate(NewID = 1:nrow(.),
         `96Well_Row` = rep(LETTERS[1:8], times = 7)[1:Nsamples],
         `96Well_Column` = gsub("Batch", "", Batch)) %>%
  left_join(old_incorrect_IL5 %>% select(SID, `Bad Relabel (Mistake)`), by = "SID") %>%
  select(SID, `Bad Relabel (Mistake)`, NewID, `96Well_Row`, `96Well_Column`)


View(for_daniel)
!any(for_daniel$SID == "105")

# 
# write_xlsx(list(RNA_Randomization = for_daniel %>% arrange(`Bad Relabel (Mistake)` == "", `Bad Relabel (Mistake)`, SID),
#                 RNA_sort_alt = for_daniel),
#            "Output/20201216_Vaping_RNA_Randomize-v2-emptyfixed.xlsx")
# 
# system("open Output/20201216_Vaping_RNA_Randomize-v2-emptyfixed.xlsx")

for_daniel %>%
  transmute(NewID, SID) %>%
  write_tsv("./Output/20201216_coreID_to_PID.txt")


for_core <-
  for_daniel %>%
  left_join(dat_qc_RNA %>% select(SID, `Conc (ng/uL)`)) %>%
  select(NewID, `96Well_Row`, `96Well_Column`, `Conc (ng/uL)`) %>%
  set_names(c("ID", "Row", "Column", "Conc (ng/ul)"))


# write_xlsx(list(SharmaVapingRNAseq = for_core),
#            path = "Output/20201217-Vaping_RNA-forcore.xlsx")

