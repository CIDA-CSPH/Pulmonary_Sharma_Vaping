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
betas <- read_tsv(here("DataProcessed/methylation/methylation_betas_BMIQ.txt"))

#Mvals
mvals <- read_tsv(here("DataProcessed/methylation/methylation_mvals_BMIQ.txt"))

#Autosomal betas
betas_auto <- read_tsv(here('DataProcessed/methylation/autosomal_betas_BMIQ.txt'))

#Metadata with autosomal intensities and predicted sex
metadata_sex <- read_csv(here("DataProcessed/methylation/metadata_all_sex_2022_08_30.csv"))

# ######################## Read in Raw IDAT Files #######################
folder_raw_dat <- searchIDATprefixes(here("DataRaw/methylation/RawIdat/"))
methylation_raw <- lapply(folder_raw_dat,readIDATpair)

#Prep Noob
methylation_noob <- openSesame(methylation_raw, prep = "QCDPB", 
                               func = NULL, BPPARAM = BiocParallel::MulticoreParam(2))

#Prep Noob + BMIQ
methylation_BMIQ <- openSesame(methylation_raw, prep = "QCDPBM", 
                               func = NULL, BPPARAM = BiocParallel::MulticoreParam(2))

#Prep Noob + BMIQ + Dye Bias (again?)
methylation_BMIQ_DB2 <- openSesame(methylation_raw, prep = "QCDPBMD", 
                                   func = NULL)

# ################## Get QC Metrics  #################
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
# write_csv(methylation_qc, here("DataProcessed/methylation/methylation_qc_metrics.csv"))

#Read in Methylation QC File
methylation_qc <- read_csv(here("DataProcessed/methylation/methylation_qc_metrics.csv"))

methylation_qc <- left_join(methylation_qc, metadata_sex %>% select(rna_id, methylation_id, sid, sentrix_name, vape_6mo_lab), by = "sentrix_name")

################## Get bad sample #################
bad_id <- metadata_sex[is.na(metadata_sex$vape_6mo_lab), "sentrix_name"] %>% as.character()

#Drop from samples
drop_badid <- function(data, badid){
  new_dat <- data[,colnames(data) != badid]
  return(new_dat)
}

#drop from betas and mvals
betas <- drop_badid(betas, bad_id)

mvals <- drop_badid(mvals, bad_id)

betas_auto <- drop_badid(betas_auto, bad_id)

m_auto <- drop_badid(m_auto, bad_id)

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

################## Dye-Bias Correction Visualization ##############
#Randomly sample three samples to check
set.seed(404)
rand_sdf <- sample(methylation_qc$sentrix_name, 3)


##############Functions needed to Calculate quantiles ############
noMasked <- function(sdf) { # filter masked probes
  sdf[!sdf$mask,,drop=FALSE]
}

InfIR <- function(sdf) {
  sdf[sdf$col == "R",,drop=FALSE]
}

InfIG <- function(sdf) {
  sdf[sdf$col == "G",,drop=FALSE]}

#Funtion to calculate what percentage of points < avg distance from line
p_inline <- function(sdf) {
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
 
 # Average distance from line
 quants_diff_mean <- mean(abs(quants_diff))

 #less than avg distance from line?
 close_to_line <- quants_diff <= quants_diff_mean & quants_diff >= -quants_diff_mean

 return(sum(close_to_line)/length(quants$x))
}
p_inline(methylation_raw[["205310760184_R08C01"]])
#Get p_inline for all samples
p_inline_res <- function(samp_names){
  deviation_df <- data.frame("Raw" = rep(NA, 3),
                           "Noob" = rep(NA, 3),
                           "Noob + BMIQ" = rep(NA, 3),
                           row.names = samp_names)
  for (h in samp_names) {
    print(h)
    raw_res <- p_inline(methylation_raw[[h]])
    noob_res <- p_inline(methylation_noob[[h]])
    bmiq_res <- p_inline(methylation_BMIQ[[h]])

    deviation_df[h,] <- c(raw_res, noob_res, bmiq_res)
  }
  
  deviation_df <- transpose(deviation_df)
  rownames(deviation_df) <- c("Raw", "Noob", "Noob + BMIQ")
  colnames(deviation_df) <- samp_names

  return(deviation_df)
}
#Get the % of samples with > avg distance from the line for each subject and each normalization procedure
deviation_res <- p_inline_res(rand_sdf) 

#write_csv(deviation_res, here("DataProcessed/methylation/dye_bias_quality_check.csv"))

# Make redgrnqq plot matrix 
plot_redgrnqq <- function(raw_dat, noob_dat, bmiq_dat, db2dat){
  
   par(mfrow=c(4,3), mar=c(3,3,2,1))
  
  #Make Raw Plots
  for (h in rand_sdf){
    sesameQC_plotRedGrnQQ(raw_dat[[h]], main= paste0(h, " (Raw)"))
    
    
  }
  #Make noob Plots
  for (i in rand_sdf){
    sesameQC_plotRedGrnQQ(noob_dat[[i]],
                              main= paste0(i, " (noob)"))
  }

  #Make noob + BMIQ Plots
  for (j in rand_sdf){
    sesameQC_plotRedGrnQQ(bmiq_dat[[j]],
                              main= paste0(j, " (noob + BMIQ)"))
  }
  #Make noob + BMIQ Plots with extra db adjustment
  for (k in rand_sdf){
    sesameQC_plotRedGrnQQ(bmiq_dat[[k]],
                          main= paste0(k, " (noob + BMIQ + db)"))
  }
}

plot_redgrnqq(raw_dat = methylation_raw,
              noob_dat = methylation_noob,
              bmiq_dat = methylation_BMIQ,
              db2dat = methylation_BMIQ_DB2)
############### Background Subtraction ########################
# par(mfrow=c(2,1), mar=c(3,3,2,1))
# sesameQC_plotBetaByDesign(methylation_raw$`205310770060_R01C01`, main="R01C01 Before", xlab="\beta")
# sesameQC_plotBetaByDesign(prepSesame(methylation_raw$`205310770060_R01C01`, prep = "QCDPB"), main="R01C01 After", xlab="\beta")

############### BMIQ + Noob ##########################
# Make redgrnqq plot matrix 
plot_betas_compare <- function(raw_dat, noob_dat, bmiq_dat, db2dat){
  
  par(mfrow=c(4,3), mar=c(3,3,2,1))
  
  #Make Raw Plots
  for (h in rand_sdf){
    sesameQC_plotBetaByDesign(raw_dat[[h]], main= paste0(h, " (Raw)"))
    
    
  }
  #Make noob Plots
  for (i in rand_sdf){
    sesameQC_plotBetaByDesign(noob_dat[[i]],
                          main= paste0(i, " (noob)"))
  }
  
  #Make noob + BMIQ Plots
  for (j in rand_sdf){
    sesameQC_plotBetaByDesign(bmiq_dat[[j]],
                          main= paste0(j, " (noob + BMIQ)"))
  }
  #Make noob + BMIQ Plots with extra db adjustment
  for (k in rand_sdf){
    sesameQC_plotBetaByDesign(bmiq_dat[[k]],
                          main= paste0(k, " (noob + BMIQ + db)"))
  }
}

plot_betas_compare(raw_dat = methylation_raw,
              noob_dat = methylation_noob,
              bmiq_dat = methylation_BMIQ,
              db2dat = methylation_BMIQ_DB2)
######################## Beta and M-Value distributions #######################

####### Betas ############
betas_long <- betas %>%
  pivot_longer(cols = !CpG_Site, names_to = c("Barcode", "Position"), values_to = "Betas", names_sep = "_") %>% 
  dplyr::mutate(Sentrix_id = paste0(Barcode, "_", Position))


#Create Plot
betas_dist <- betas_long %>% 
  ggplot(aes(x = Betas, group = Sentrix_id, col = Barcode)) +
  geom_density(aes(y = ..scaled..)) +
  labs(x = expression(beta),
       y = "Density")+
  theme(legend.position = "bottom")


####### M-Values ############
mvals_long <- mvals %>%
  pivot_longer(cols = !CpG_Site, names_to = c("Barcode", "Position"), values_to = "M", names_sep = "_") %>% 
  dplyr::mutate(Sentrix_id = paste0(Barcode, "_", Position))

#Create Plot
mvals_dist <- mvals_long %>% 
  ggplot(aes(x = M, group = Sentrix_id, col = Barcode)) +
  geom_density(aes(y = ..scaled..)) +
  labs(x = "M",
       y = "Density")+
  theme(legend.position = "bottom",
        axis.title.y = element_blank())


########### Arrange Plots ##########
ggarrange(betas_dist, mvals_dist, nrow = 1, ncol = 2, common.legend = T, legend = "bottom")

######################### Sex Check ##################################
#Clinical Sex
clin_sex_plot <- metadata_sex %>% 
  ggplot(aes(x = log2(medianX), y = log2(medianY), col = sex_lab)) +
  geom_text(aes(label = sid)) +
  labs(col = "Sex") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("deepskyblue", "deeppink3"))

#Predicted Sex
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


make_mds_plots(betas, metadata_sex, plot = "vape")
make_mds_plots(betas_auto, metadata_sex)
