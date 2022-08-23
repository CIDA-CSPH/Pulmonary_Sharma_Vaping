####################Load Libraries########################
library(tidyverse)
library(here)
library(kableExtra)
library(sesame)
library(readxl)
library(parallel)
library(randomForest)
library(janitor)
library(minfi)
library(data.table)
library(ggpubr)

kablize <- function(tab, digits = 3) {
  kable(tab,digits = digits, booktabs = T) %>% 
    kableExtra::kable_styling(latex_options = c("scale_down", "striped"), position = "center")
}

format_num <- function(number, digits = 0, format = "f") {
  formatC(number, digits = digits, format = format, big.mark = ",")
}

################# Read in Betas, M-Values, and Metadata ##################
#betas
betas <- fread(here("DataProcessed/methylation/methylation_betas.txt"))

#Mvals
mvals <- fread(here("DataProcessed/methylation/methylation_mvals.txt"))

#Metadata with autosomal intensities and predicted sex
metadata_sex <- read_csv(here("DataProcessed/methylation/metadata_all_sex_2022_08_19.csv")) %>% 
  dplyr::mutate(sex_lab = if_else(sex_lab == "Male", "Female", "Male"),
                Plate = str_split(sentrix_name, "_") %>% map_chr(.,1),
                Position = str_split(sentrix_name, "_") %>% map_chr(.,2))

######################## Beta and M-Value distributions #######################

####### Betas ############
# betas_long <- betas %>%
#   pivot_longer(cols = !CpG_Site, names_to = "Sentrix_id", values_to = "Betas" ) %>%
#   mutate(Plate = str_split(Sentrix_id, "_") %>% map_chr(.,1),
#          Position = str_split(Sentrix_id, "_") %>% map_chr(.,2))

betas_long <- fread(here("DataProcessed/methylation/methylation_betas_long.txt"))

#Create Plot
betas_dist <- betas_long %>% 
  ggplot(aes(x = Betas, group = Sentrix_id, col = Plate)) +
  geom_density(aes(y = ..scaled..)) +
  labs(x = expression(beta),
       y = "Density")+
  theme(legend.position = "bottom")


####### M-Values ############
# mvals_long <- betas %>%
#   pivot_longer(cols = !CpG_Site, names_to = "Sentrix_id", values_to = "M" ) %>%
#   mutate(Plate = str_split(Sentrix_id, "_") %>% map_chr(.,1),
#          Position = str_split(Sentrix_id, "_") %>% map_chr(.,2))

m_long <- fread(here("DataProcessed/methylation/methylation_mvals_long.txt"))

#Create Plot
mvals_dist <- m_long %>% 
  ggplot(aes(x = M, group = Sentrix_id, col = Plate)) +
  geom_density(aes(y = ..scaled..)) +
  labs(x = "M",
       y = "Density")+
  theme(legend.position = "bottom",
        axis.title.y = element_blank())

########### Arrange Plots ##########
ggarrange(betas_dist, mvals_dist, nrow = 1, ncol = 2, common.legend = T, legend = "bottom")

######################### Clustering by Sex ##################################
#Clinical Sex
clin_sex_plot <- metadata_sex %>% 
  ggplot(aes(x = log2(medianX), y = log2(medianY), col = sex_lab)) +
  geom_text(aes(label = new_id)) +
  labs(col = "Clinical Sex") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("deepskyblue", "deeppink3"))

#Predicted Sex
predicted_sex_plot <- metadata_sex %>% 
  ggplot(aes(x = log2(medianX), y = log2(medianY), col = pred_sex)) +
  geom_text(aes(label = new_id)) +
  labs(col = "Predicted Sex") +
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("deepskyblue", "deeppink3"))

ggarrange(clin_sex_plot, predicted_sex_plot, nrow = 1, ncol = 2, common.legend = T)

######################### Clustering by Plate ##################################
betas <- as.data.frame(betas)

rownames(betas) <- as.list(betas$CpG_Sites)

betas_mat <- betas %>% 
  select(-CpG_Sites) %>% as.matrix()

mdsPlot(betas_mat,
        sampNames = colnames(betas_mat),
        sampGroups = metadata_sex$)

metadata_unjoined %>% 
  group_by(sex_lab) %>% 
  summarise(n = n())
