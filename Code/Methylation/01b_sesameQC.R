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

################# Read in Betas, M-Values, and Metadata ##################
#betas
betas <- read_tsv(here("DataProcessed/methylation/methylation_betas_final_2022_09_27.txt"))

#Autosomal betas
betas_auto <- read_tsv(here('DataProcessed/methylation/autosomal_betas_final_2022_09_30.txt'))

#Metadata with autosomal intensities and predicted sex
metadata_sex <- read_csv(here("DataProcessed/methylation/metadata_all_sex_2022_08_30.csv"))

# ######################## Read in Raw IDAT Files #######################
set.seed(404)
folder_raw_dat <- sample(searchIDATprefixes(here("DataRaw/methylation/RawIdat/")), 3)

methylation_raw <- lapply(folder_raw_dat,readIDATpair)

#Prep Noob + BMIQ
methylation_BMIQ <- openSesame(methylation_raw, prep = "QCDPBM", 
                               func = NULL, BPPARAM = BiocParallel::MulticoreParam(2))

#Prep Noob + BMIQ + Dye Bias (last)
methylation_BMIQ_DB2 <- openSesame(methylation_raw, prep = "QCPBMD", 
                                   func = NULL, BPPARAM = BiocParallel::MulticoreParam(2))


# Identify Outlier Sample By Mean Intensity -------------------------------
# methylation_qc <- do.call(rbind, lapply(methylation_raw, function(x)
#                             as.data.frame(sesameQC_calcStats(x, "intensity")))) %>% 
#   mutate(sentrix_name = rownames(.))
# 
# rownames(methylation_qc) <- NULL
# 
# ################## Find Outliers  #################
# outlier_quant <- quantile(methylation_qc$mean_intensity, probs = c(0.01, 0.99))
# 
# methylation_qc <- methylation_qc %>%
# mutate(outlier = if_else(mean_intensity < outlier_quant[1], T, F))
# 
# #Write out methylation qc file
# write_csv(methylation_qc, here("DataProcessed/methylation/QC/methylation_qc_metrics.csv"))

#Read in Methylation QC File
methylation_qc <- read_csv(here("DataProcessed/methylation/QC/methylation_qc_metrics.csv"))

methylation_qc <- left_join(methylation_qc, metadata_sex %>% select(rna_id, methylation_id, sid, sentrix_name, vape_6mo_lab), by = "sentrix_name")


# Remove Sample Without Vape Status ---------------------------------------
bad_id <- metadata_sex[is.na(metadata_sex$vape_6mo_lab), "sentrix_name"] %>% as.character()

#Drop from samples
drop_badid <- function(data, badid){
  new_dat <- data[,colnames(data) != badid]
  return(new_dat)
}

#drop from betas and mvals
betas <- drop_badid(betas, bad_id)

betas_auto <- drop_badid(betas_auto, bad_id)

#drop from metadata and methylation_qc
methylation_qc <- methylation_qc %>% drop_na(vape_6mo_lab)

metadata_sex <- metadata_sex %>% drop_na(vape_6mo_lab)
################## Plot Outliers  #################

methylation_qc %>% 
  ggplot(aes(y = log2(mean_intensity))) +
  geom_boxplot() +
  geom_text(data = methylation_qc[methylation_qc$outlier == T, ], aes(y = log2(mean_intensity), x = -0.02, label = methylation_qc[methylation_qc$outlier == T, ]$sid), col = "red")+
  geom_hline(aes(yintercept = log2(quantile(methylation_qc$mean_intensity, probs = 0.01)), col = "1st Percentile"), linetype = "dashed") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom") +
  scale_color_manual(values = c("red")) +
  labs(y = "log2(Mean Intensity)",
       color = NULL)


# Dye Bias QC -------------------------------------------------------------
rand_sdf <- names(methylation_raw)

## Functions to get red/grn probes only --------
noMasked <- function(sdf) { # filter masked probes
  sdf[!sdf$mask,,drop=FALSE]
}

InfIR <- function(sdf) {
  sdf[sdf$col == "R",,drop=FALSE]
}

InfIG <- function(sdf) {
  sdf[sdf$col == "G",,drop=FALSE]}


# Get the mean distance from the x=y for unnormalized data to compare --------
#All Green Probes
dG <- InfIG(noMasked(methylation_raw[[1]])) 
#All Red Probes
dR <- InfIR(noMasked(methylation_raw[[1]]))
# all red probes intensity
all_red <- c(dR$MR, dR$UR)
# all green probes intensity
all_grn <- c(dG$MG, dG$UG)
# Quantiles for each point
quants <- qqplot(all_red, all_grn, plot.it = F)

quants_diff <- quants$y - quants$x

raw_quants_mean <- abs(mean(quants_diff))

## Funtion to calculate what percentage of points < avg distance from line
p_inline <- function(sdf, mean_val) {
  #All Green Probes
  dG <- InfIG(noMasked(sdf)) 
  #All Red Probes
  dR <- InfIR(noMasked(sdf))
  # all red probes intensity
  all_red <- c(dR$MR, dR$UR)
  # all green probes intensity
  all_grn <- c(dG$MG, dG$UG)
 # Quantiles for each point
 quants <- qqplot(all_red, all_grn, plot.it = F)
 
 quants_diff <- quants$y - quants$x

 #less than avg distance from line?
 close_to_line <- quants_diff <= mean_val & quants_diff >= -mean_val

 return(sum(close_to_line)/length(quants$x))
}

#Get p_inline for all samples
p_inline_res <- function(samp_names, raw_mean){
  deviation_df <- data.frame("Raw" = rep(NA, 3),
                           "dbNL + Noob + BMIQ" = rep(NA, 3),
                           "Noob + BMIQ + dbNL" = rep(NA, 3),
                           row.names = samp_names)
  for (h in samp_names) {
    print(h)
    raw_res <- p_inline(methylation_raw[[h]], mean_val = raw_mean)
    bmiq_res <- p_inline(methylation_BMIQ[[h]], mean_val = raw_mean)
    bmiq_db2_res <- p_inline(methylation_BMIQ_DB2[[h]], mean_val = raw_mean)

    deviation_df[h,] <- c(raw_res, bmiq_res, bmiq_db2_res)
  }
  
  deviation_df <- transpose(deviation_df)
  rownames(deviation_df) <- c("Raw", "dbNL + Noob + BMIQ", "Noob + BMIQ + dbNL")
  colnames(deviation_df) <- samp_names

  return(deviation_df)
}
#Get the % of samples with > avg distance from the line for each subject and each normalization procedure
deviation_res <- p_inline_res(rand_sdf, raw_quants_mean) 

#write_csv(deviation_res, here("DataProcessed/methylation/QC/dye_bias_quality_check.csv"))

# Make redgrnqq plot matrix 
plot_redgrnqq <- function(raw_dat, bmiq_dat, db2dat){
  
   par(mfrow=c(3,3), mar=c(4,4,2,1))
  
  #Make Raw Plots
  for (h in rand_sdf){
    sesameQC_plotRedGrnQQ(raw_dat[[h]], main= paste0(h, " (Raw)"))
    
    
  }
  #Make dbNL + noob + BMIQ Plots
  for (i in rand_sdf){
    sesameQC_plotRedGrnQQ(bmiq_dat[[i]],main= paste0(i, " (dbNL + noob + BMIQ)"))
  }

  #Make noob + BMIQ +dbNL Plots
  for (j in rand_sdf){
    sesameQC_plotRedGrnQQ(db2dat[[j]],main= paste0(j, " (noob + BMIQ + dbNL)"))
  }
}

plot_redgrnqq(raw_dat = methylation_raw,
              bmiq_dat = methylation_BMIQ,
              db2dat = methylation_BMIQ_DB2)


# Check Betas Dist'n ------------------------------------------------------
## This function makes a matrix of plots for differnet normalizations
plot_betas_compare <- function(raw_dat, bmiq_dat, db2dat){
  
  par(mfrow=c(3,3), mar=c(3,3,2,1))
  
  #Make Raw Plots
  for (h in rand_sdf){
    sesameQC_plotBetaByDesign(raw_dat[[h]], main= paste0(h, " (Raw)"))
    
    
  }
  #Make dbNL + noob + BMIQ Plots
  for (i in rand_sdf){
    sesameQC_plotBetaByDesign(bmiq_dat[[i]],
                          main= paste0(i, " (dbNL + noob + BMIQ)"))
  }
  
  #Make noob + BMIQ +dbNL Plots
  for (j in rand_sdf){
    sesameQC_plotBetaByDesign(db2dat[[j]],
                          main= paste0(j, " (noob + BMIQ + dbNL)"))
  }
}

plot_betas_compare(raw_dat = methylation_raw,
              bmiq_dat = methylation_BMIQ,
              db2dat = methylation_BMIQ_DB2)


# Beta and M-Value Dist'ns (all samples) -------------------------------------
## Get long-format betas
betas_long <- betas %>%
  pivot_longer(cols = !CpG_Site, names_to = c("Barcode", "Position"), values_to = "Betas", names_sep = "_") %>% 
  dplyr::mutate(Sentrix_id = paste0(Barcode, "_", Position))


## Create Dist'n plot
betas_dist <- betas_long %>% 
  ggplot(aes(x = Betas, group = Sentrix_id, col = Barcode)) +
  geom_density(aes(y = ..scaled..)) +
  labs(x = expression(beta),
       y = "Density")+
  theme(legend.position = "bottom")


## Get long format mvals
mvals <- BetaValueToMValue(betas)

mvals_long <- mvals %>%
  pivot_longer(cols = !CpG_Site, names_to = c("Barcode", "Position"), values_to = "M", names_sep = "_") %>% 
  dplyr::mutate(Sentrix_id = paste0(Barcode, "_", Position))

## Create Dist'n plot
mvals_dist <- mvals_long %>% 
  ggplot(aes(x = M, group = Sentrix_id, col = Barcode)) +
  geom_density(aes(y = ..scaled..)) +
  labs(x = "M",
       y = "Density")+
  theme(legend.position = "bottom",
        axis.title.y = element_blank())


## Create single visual
ggarrange(betas_dist, mvals_dist, nrow = 1, ncol = 2, common.legend = T, legend = "bottom")


# Visualization for checking sex ------------------------------------------

## Clinical Sex
clin_sex_plot <- metadata_sex %>% 
  ggplot(aes(x = log2(medianX), y = log2(medianY), col = sex_lab)) +
  geom_text(aes(label = sid)) +
  labs(col = "Sex") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("deepskyblue", "deeppink3"))

## Predicted Sex
predicted_sex_plot <- metadata_sex %>% 
  ggplot(aes(x = log2(medianX), y = log2(medianY), col = pred_sex)) +
  geom_text(aes(label = sid)) +
  labs(col = "Predicted Sex") +
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("deepskyblue", "deeppink3"))

ggarrange(clin_sex_plot, predicted_sex_plot, nrow = 1, ncol = 2, common.legend = T)

######################### Whole Genome mdsPlots ##################################

#convert cpgsites back to rownames
mds_prep <- function(dat) {
  dat <- as.data.frame(dat)
  rownames(dat) <- dat[,"CpG_Site"]
  new_dat <- dat %>% 
    select(-CpG_Site) %>% 
    as.matrix()
  return(new_dat)
}

#Code to prep data and make mds plots
make_mds_plots <- function(vals, meta, plot = c("sex", "vape", "center")){
  #Get only samples with methylation data and arrange in same order as betas so that you can use for mds plots
  meta <- meta %>% 
    filter(sentrix_name %in% colnames(vals)) %>% 
    arrange(factor(sentrix_name, levels = colnames(vals)))
  
  #Prep values for mds plots
  vals <- mds_prep(vals)
  
  #Get mvals
  mvals <- BetaValueToMValue(vals)
  
  # Set Color Pallattes
  sex_pal <- c("deepskyblue", "deeppink3")
  center_colors <- c("cyan3", "chartreuse3", "darkorange")
  
  #Number of Sites to use
  num_1000 = 1000
  num_10000 = 10000
  
  
  #MDSPlots (Betas)
  par(mfrow = c(1,2))
  if ("sex" %in% plot) {
  #Clustering by Sex
  mdsPlot(vals, sampGroups = meta$sex_lab, sampNames = meta$sid,
          numPositions=num_1000, ylim = c(-10,10), main = "Sex - Betas", pal = sex_pal)
  
  mdsPlot(mvals, sampGroups = meta$sex_lab, sampNames = meta$sid,
          numPositions=num_1000, ylim = c(-10,10), main = "Sex - M-values", pal = sex_pal)
  }
  
  if ("vape" %in% plot) {
  #Clustering by Vape Status
  mdsPlot(vals, sampGroups = meta$vape_6mo_lab, sampNames = meta$sid,
          numPositions=num_1000, ylim = c(-10,10), main = "Vape Status - Betas")
  
  mdsPlot(mvals, sampGroups = meta$vape_6mo_lab, sampNames = meta$sid,
          numPositions=num_1000, ylim = c(-10,10), main = "Vape Status - M-values")
  }
  
  if ("center" %in% plot) {
  #Clustering by center
  mdsPlot(vals, sampGroups = meta$recruitment_center, sampNames = meta$sid,
          numPositions=num_1000, ylim = c(-10,10), pal = center_colors, main = "Recruitment Center - Betas")
  
  mdsPlot(mvals, sampGroups = meta$recruitment_center, sampNames = meta$sid,
          numPositions=num_1000, ylim = c(-10,10), pal = center_colors, main = "Recruitment Center - M-values")
  }
}

