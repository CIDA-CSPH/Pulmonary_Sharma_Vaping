#load in libraries
library(readr)
library(arsenal)
library(tidyverse)
library(janitor)
library(here)
library(forcats)

#read in data
all_metadata <- readr::read_csv(here::here("choo_code_archive/AdolescentVapingRNASeq/Data/20210325/all_data_Latino_Youth_survey_clean_pulmfxn_FINAL_3.23.21.csv"))

head(all_metadata)
glimpse(all_metadata)

#clean variable names
vape_dat_clean <- clean_names(all_metadata)

#Fix Age
vape_dat_clean <- vape_dat_clean %>% 
  mutate(age = age + 7)

#Check cities
vape_dat_clean %>% 
  group_by(city) %>% 
  summarise(N = n()) %>% 
  ungroup

#Fix cities
vape_dat_clean <- vape_dat_clean %>% 
mutate(recruitment_center = case_when(city == "Aurora" ~ "Aurora",
                                      city == "Pueblo" | city == "Avondale" | city == "Mineral" ~ "Pueblo",
                                      T ~ "CommCity/Denver"))

#Check 'sex' (gender variable is actually sex variable)
vape_dat_clean %>% 
  group_by(gender) %>%
  summarise(N = n())

vape_dat_clean %>% 
  filter(gender==5) %>% 
  select(gender,gender_txt)

#Fix 'Gender' -> 'Sex'
vape_dat_clean <- vape_dat_clean %>% 
  mutate(sex_lab = case_when(gender == 1 ~ 'Male',
                             gender == 2 ~ 'Female',
                             gender == 5 ~ 'Female')) #Cheyret used Methylation data to confirm non-binary indiv is Female
#Check Grade
vape_dat_clean %>% 
  group_by(grade) %>% 
  summarise(N = n())

#Add Grade Label column
vape_dat_clean <- vape_dat_clean %>% 
  mutate(grade_lab = factor(grade, labels = c('7th', '8th','Freshman','Sophomore','Junior','Senior'), ordered = is.ordered(grade) ))

vape_dat_clean$grade_lab

#check Latino
vape_dat_clean %>% 
  group_by(latino) %>% 
  summarise(N = n())

#Latino Label
vape_dat_clean <- vape_dat_clean %>% 
  mutate(latino_lab = if_else(latino == 1, 'LatinX', 'Non-LatinX'))

#check ever_vaped?
vape_dat_clean %>% 
  group_by(ever_vape) %>% 
  summarise(N = n()) %>% 
  mutate(Description = if_else(ever_vape == 0, 'never vaped', "vaped")) %>% 
  ungroup()

#Make ever_vaped column
vape_dat_clean <- vape_dat_clean %>% 
  mutate(ever_vape_lab = if_else(ever_vape == 0, "Never Vaped", "Vaped"))

#Check vaped in last 30 days?
vape_dat_clean %>% 
  mutate(vape_days = if_else(ever_vape == 0 & is.na(vape_days), 0, vape_days)) %>% 
  group_by(vape_days > 0) %>% 
  summarise(N = n())

#Fix NA's for vaped in 30 days
vape_dat_clean <- vape_dat_clean %>% 
  mutate(vape_days = if_else(ever_vape == 0 & is.na(vape_days), 0, vape_days))

#Make Vaped in 30 days label column
vape_dat_clean <- vape_dat_clean %>% 
  mutate(vape_30_lab = if_else(vape_days > 0, 'Vaped in last 30 Days', 'Did Not Vape in Last 30 Days'))

#Check Vaped in 6mo?
vape_dat_clean %>% 
  group_by(last_vape < 5 | last_vape > 8) %>% 
  summarise(N = n()) %>% 
  ungroup()

#create column fixing NA's
vape_dat_clean <- vape_dat_clean %>% 
  mutate(vape_6mo = case_when(last_vape < 5 | last_vape > 8 ~ TRUE,
                              last_vape %in% c(5,6,7,8) ~ FALSE,
                              ever_vape == 0 & is.na(last_vape) ~ FALSE,
                              vape_days > 0 ~ TRUE))

#Create labeled column for vaped in 6 mo
vape_dat_clean <- vape_dat_clean %>% 
  mutate(vape_6mo_lab = if_else(vape_6mo == TRUE, 'Vaped in Last 6 Months', 'Did Not Vape in Last 6 Months'))

#subset data for table1
tab1_dat <- vape_dat_clean %>% 
  select(sid, recruitment_center, age, sex_lab, grade_lab, latino_lab, ever_vape_lab, vape_30_lab, vape_6mo_lab, fev1, fvc, r5, x20, fev1_fvc)

#output csv
# write_csv(tab1_dat, here('DataProcessed/metadata_cleaning/table1_clean_data_yyyy_mm_dd.csv'))


