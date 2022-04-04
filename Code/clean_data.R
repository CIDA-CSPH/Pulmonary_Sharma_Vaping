df.all <- read.csv("Without_new_variables_Latinx_Youth_Vaping_Project_8.20.2020.csv")
df.all2 <- read.csv("new_variables_Latinx_Youth_Vaping_Project_8.20.2020.csv")
library(dummies)
library(splitstackshape)
library(dplyr)

head(df$Ejuice_strength1)
df <- df.all %>%
  select(SID, Race, ever_tobacco, device_type, Ejuice_strength1, ejuice_flavor, Household_use, 
         friend1_typeact, friend2_typeact, friend3_typeact, friend4_typeact, friend5_typeact, friend6_typeact)

vars_to_split <- c("Race", "ever_tobacco", "device_type", "Ejuice_strength1", "ejuice_flavor", "Household_use", 
                   "friend1_typeact", "friend2_typeact", "friend3_typeact", "friend4_typeact", "friend5_typeact", 
                   "friend6_typeact")
table(df$ever_tobacco)

for(i in vars_to_split){
  df <- cSplit_e(df, i, sep = ",", fill = 0)
}


df <- df %>%
  select(-Race_6, -ever_tobacco_3)

new.names <- c(colnames(df)[1:13], 
                  "Race_AIAN", "Race_Asian", "Race_black", "Race_NH_PI", "Race_white", "Race_other",
                  "Ever_cig", "Enver_cigars", "Ever_smokeless", "Ever_Snus", "Ever_hookah", "Ever_bidis", "Ever_vape",
                  "device_JUUL","device_Njoy","device_pen","device_ecig","device_mods",
                  "device_hookahpen","device_ehookah","device_ecigar","device_epipe","device_other",
                  colnames(df)[37:50],
               "flavor_tob", "flavor_mental", "flavor_clove", "flavor_fruit", "flavor_choc", "flavor_alcoholic",
               "flavor_nonalc", "flavor_candy", "flavor_other",
               "Home_cig", "Home_cigar", "Home_chew", "Home_hookah", "Home_ecig", "Home_noone",
               "friend1_activity", "friend1_hangout", "friend1_close", "friend1_alike", "friend1_other",
               "friend2_activity", "friend2_hangout", "friend2_close", "friend2_alike", "friend2_other",
               "friend3_activity", "friend3_hangout", "friend3_close", "friend3_alike", "friend3_other",
               "friend4_activity", "friend4_hangout", "friend4_close", "friend4_alike", "friend4_other",
               "friend5_activity", "friend5_hangout", "friend5_close", "friend5_alike", "friend5_other",
               "friend6_activity", "friend6_hangout", "friend6_close", "friend6_alike", "friend6_other")

new.names <- data.frame(old_name = colnames(df), 
                        new_name = new.names)

colnames(df) <- new.names$new_name

df <- df %>%
  select(-all_of(vars_to_split))

new.df.all <- inner_join(df.all, df, by = "SID")
write.csv(new.df.all, "all_data.csv")
