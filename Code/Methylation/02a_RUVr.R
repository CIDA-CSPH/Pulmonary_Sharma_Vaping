# Load Libraries ----------------------------------------------------------
library(tidyverse)
library(sesame)
library(sesameData)
library(minfi)
library(missMethyl)
library(here)
library(limma)
library(ggrepel)
library(ggthemes)
library(RUVSeq)
library(RColorBrewer)

# Load in Required Data ---------------------------------------------------
## targets to read in ---------------------------------------------------
targets <- read_csv(here("DataProcessed/methylation/minfi_targets.csv"))

## other clinical metadata ---------------------------------------------------
clin_metadata <- read_csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv"))

clin_metadata <- clin_metadata %>% 
  drop_na(methylation_id, vape_6mo_lab)

## M Values ----------------------------------------------------------------
mvals <- read_tsv(here("DataProcessed/methylation/methylation_mvals_final_2022_09_27.txt")) %>% as.data.frame()

rownames(mvals) <- mvals$CpG_Site

## Drop the NA Vaper
mvals <- mvals %>% select(-CpG_Site, -targets$sentrix_name[is.na(targets$vape_6mo_lab)])

# Drop the NA Vaper From Targets-------------------------------------------------------
targets <- targets %>% drop_na(vape_6mo_lab)


# Elbow method for finding optimal ruv clusters ---------------------------
set.seed(404)

#get scaled data
scaled_data <- as.matrix(scale(mvals))

# Compute and plot wss for k = 2 to k = 5.
k.max <- 5
wss <- sapply(1:k.max, 
              function(k){kmeans(scaled_data, k, nstart = 10, iter.max = 15)$tot.withinss})

elbow_tib <- tibble(k = seq(1,5,1), 
                    wss = wss)


elbow_plot <- elbow_tib %>% 
  ggplot(aes(x = k, y = wss)) +
  geom_point() +
  geom_line() +
  labs(x = "K", y = "Total Within-Clusters Sum of Squares", title = "Elbow Plot", col = "") + 
  scale_x_continuous(breaks = seq(1,10,1))

elbow_plot


clin_metadata <- clin_metadata[match(names(mvals), clin_metadata$sentrix_name),]
# #Continue with K = 2
# plan(multisession)
# # Get the first pass residuals
# # preruv_residuals <- 
# #   future_apply(mvals, 1,
# #         function(y) { lm(y ~ vape_6mo_lab + sex_lab + age,
# #                          clin_metadata) %>% 
# #             residuals()}) %>% 
#   t
# 
# colnames(preruv_residuals) <- names(mvals)

#preruv_residuals_write <- preruv_residuals %>% as.data.frame %>% rownames_to_column(var = "CpG_Site")
#write_csv(preruv_residuals_write, here("DataProcessed/methylation/ruv_residuals_2022_09_27.csv"))

preruv_residuals <- read_csv(here("DataProcessed/methylation/ruv_residuals_2022_09_27.csv")) %>% as.data.frame()

rownames(preruv_residuals) <- preruv_residuals$CpG_Site

preruv_residuals <- preruv_residuals %>% select(-CpG_Site)

get_ruv_res <- function (vals, resid, k) {
  #Run RUVr
  ruv_res <- RUVr(
    x = vals %>% as.matrix(),
    k = k,
    residuals = resid %>% as.matrix(),
    isLog = T)
  
  #RUN PCA
  pca_results <-
    prcomp(ruv_res$normalizedCounts %>% t, center = T, scale = F)
  
  #Get the percent explained
  pca_percent_exp <-
    (pca_results$sdev^2) %>% `/`(., sum(.)) %>% `*`(100) %>% formatC(digits = 2) %>%
    gsub(" ", "", .) %>%
    paste0(" (", ., "%)") %>%
    paste0("PC", 1:length(.), .)
  
  #Bind the datframe for plotting
  pca_df <-
    clin_metadata %>%
    bind_cols(pca_results$x[ , 1:4] %>% as_tibble)
  
  return(pca_df)
}


ruv_k0 <- get_ruv_res(mvals, preruv_residuals, k = 0)

ruv_k2 <- get_ruv_res(mvals, preruv_residuals, k = 2)



#Make Plots

vape_cols <- brewer.pal(3, "Dark2")

pca_plots <- function(ruv_res, plot) {

  if (plot == "vape")  {
    # Vape Status
    return(ruv_res %>% 
      ggplot(aes(PC1, PC2, color = vape_6mo_lab)) + #, color = sex, shape = sex)) +
      geom_point(alpha = 0.75) +
      geom_text_repel(
        aes(label = sid),
        color = "black", size = 2) +
      #xlab(pca_percent_exp[1]) + ylab(pca_percent_exp[2]) +
      ggtitle("Methylation PCA (Vape)") +
      theme_few() +
      geom_blank()+
      scale_color_manual(name = "Vape Status", labels = c("No Vape", "Vape"), values = vape_cols[1:2]))
  }
  if (plot == "center") {
    # Recruitment Center
    return(ruv_res %>% 
      ggplot(aes(PC1, PC2, color = recruitment_center)) +
      geom_point(alpha = 0.75) +
      geom_text_repel(
        aes(label = sid),
        color = "black", size = 2) +
      #xlab(pca_percent_exp[1]) + ylab(pca_percent_exp[2]) +
      ggtitle("Methylation PCA (Center)") +
      theme_few() +
      geom_blank())
  }
  if (plot == "sex") {
    #Sex
    return(ruv_res %>% 
      ggplot(aes(PC1, PC2, color = sex_lab)) +
      geom_point(alpha = 0.75) +
      geom_text_repel(
        aes(label = sid),
        color = "black", size = 2) +
      #xlab(pca_percent_exp[1]) + ylab(pca_percent_exp[2]) +
      ggtitle("Methylation PCA (Sex)") +
      theme_few() +
      geom_blank())
  }
}

pca_plots(ruv_k0, plot = "vape")
pca_plots(ruv_k2, plot = "vape")

pca_plots(ruv_k0, plot = "center")
pca_plots(ruv_k2, plot = "center")

pca_plots(ruv_k0, plot = "sex")
pca_plots(ruv_k2, plot = "sex")
