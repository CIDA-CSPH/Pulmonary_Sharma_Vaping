# Read Libraries ----------------------------------------------------------
library(tidyverse)
library(here)
library(minfi)
library(sesame)
library(skewr)
library(readxl)
library(bumphunter)
library(RColorBrewer)


# Set Base Directory ------------------------------------------------------
baseDir <- here("DataRaw/methylation/RawIdat")

# Read in Metadata and Comparison -----------------------------------------------------
## clinical metadata
clin_metadata <- read_csv(here("DataProcessed/clinical_metadata/master_clinical_metadata_2022_09_02.csv"))

#betas
nb_betas <- read_tsv(here("DataProcessed/methylation/methylation_betas_BMIQ.txt"))

#Mvals
nb_mvals <- read_tsv(here("DataProcessed/methylation/methylation_mvals_BMIQ.txt"))

#Autosomal betas
nb_betas_auto <- read_tsv(here('DataProcessed/methylation/autosomal_betas_BMIQ.txt'))

#clinical metadata with chromosomal intensity values
metadata_sex <- read_csv(here("DataProcessed/methylation/metadata_all_sex_2022_08_30.csv"))
#Label column for plots 
metadata_sex <- metadata_sex %>% 
  mutate(vape_text = if_else(vape_6mo_lab == "Vaped in Last 6 Months", "Vaped", "Not Vaped"))
#Read in Methylation QC File
methylation_qc <- read_csv(here("DataProcessed/methylation/methylation_qc_metrics.csv"))

methylation_qc <- left_join(methylation_qc, metadata_sex %>% select(rna_id, methylation_id, sid, sentrix_name, vape_6mo_lab), by = "sentrix_name")

# Read in Sample Data -----------------------------------------------------

#Read in and join methyl metadata
samplesheet <- left_join(
  read_csv(here("DataRaw/methylation/Metadata/SSharma_48_samplesheet_05202021.csv"), skip = 7) %>%
    select(-`Sample_Group`, -`Pool_ID`),
  read_excel(here("DataRaw/methylation/Metadata/Sharma_48_5192021Controls.xlsx")) %>%
    mutate(`Sample Name` = as.numeric(`Sample Name`)),
  by = c("Sample_Name" = "Sample Name"))

samplesheet <- samplesheet %>% 
  mutate(sample_join = paste0("Sample", str_pad(as.character(Sample_Name), 2, pad = "0")),
         sentrix_name = paste0(`Sentrix Barcode`, "_", `Sentrix Position`))

# Getting Filename prefixes for readIDATpairs
folder_raw_dat <- here("DataRaw/methylation/RawIdat/")

targets <- tibble(file_path =
                    list.files(path = folder_raw_dat, pattern = "*.idat|*.IDAT",
                               full.names = F, recursive = T)) %>%
  mutate(`Sentrix Barcode` = str_split(file_path, pattern = "/|_") %>% map_chr( ~ .[4]), # chip
         `Sentrix Position` = str_split(file_path, pattern = "/|_") %>% map_chr( ~ .[5])) %>% # row
  left_join(samplesheet,
            by = c("Sentrix Barcode", "Sentrix Position")) %>%
  mutate(file_path_prefix = gsub("_Red.idat|_Grn.idat", "", file_path) %>%
           paste0(folder_raw_dat, .)) %>%
  filter(!duplicated(file_path_prefix))

targets <- targets %>% 
  dplyr::rename(Basename = file_path_prefix)

#join metadata back to targets
targets <- left_join(targets, clin_metadata %>% 
                       select(methylation_id,sid,sex_lab, age, recruitment_center, vape_6mo_lab),
                     by = c("sample_join" = "methylation_id"))


# Read in Methylation Data -----------------------------------------------------
RGset <- read.metharray.exp(targets = targets)

getManifest(RGset)

# Subset Probes to be the same as Noob + BMIQ -----------------------------
RGset_compare <- subsetByLoci(RGset, includeLoci = nb_betas$CpG_Site, keepSnps = F)

# Perform normalization
mset_raw <- preprocessRaw(RGset_compare)
mset_noob <- preprocessNoob(RGset_compare, dyeCorr = T, dyeMethod = "reference")
mset_noob_swan <- preprocessSWAN(RGset_compare, mSet = mset_noob)
ns_betas <- getBeta(mset_noob_swan)

# Dye Bias  ---------------------------------------------------------------


# Beta Dist'n by Probe Type -----------------------------------------------
set.seed(404)
rand_samp <- sample(methylation_qc$sentrix_name, 3)

par(mfrow = c(2,3))
for (i in rand_samp) {
plotBetasByType(mset_noob[,i], 
                main = paste0(i, " Noob"))
}
for (j in rand_samp) {
  plotBetasByType(mset_noob_swan[,j],
                  main = paste0(j, " Noob + SWAN"))
}


# Beta Dist'n by Probe Type and Color -------------------------------------

## Source Code ------------------------------------
plotBetasByTypeColor <- function(data, probeTypes = NULL, legendPos = "top",
                            colors = c("black", "red", "green", "blue"), main = "",
                            lwd = 3, cex.legend = 1) {
  # Check inputs
  if (!(is(data, "MethylSet") || is.matrix(data) || is.vector(data) ||
        is(data, "DelayedMatrix"))) {
    stop("'data' needs to be a 'MethylSet' or a matrix, 'DelayedMatrix', ",
         "or vector of beta values")
  }
  if (!is.vector(data)) {
    if (ncol(data) > 1) stop("'data' must only contain one sample")
  }
  if (!is(data, "MethylSet")) {
    r <- range(data, na.rm = TRUE)
    if (!(r[1] >= 0 && r[2] <= 1)) {
      stop("'data' should be beta values in the range [0, 1]")
    }
  }
  if (is(data, "MethylSet")) {
    if (is.null(probeTypes)) {
      typeI_grn <- getProbeInfo(data, type = "I-Green")[, c("Name", "nCpG")]
      typeI_red <- getProbeInfo(data, type = "I-Red")[, c("Name", "nCpG")]
      typeII <- getProbeInfo(data, type = "II")[, c("Name", "nCpG")]
      probeTypes <- rbind(typeI_grn, typeI_red, typeII)
      probeTypes$Type <- rep(
        x = c("I-grn","I-red", "II"),
        times = c(nrow(typeI_grn), nrow(typeI_red), nrow(typeII)))
    }
  }
  if (!all(c("Name", "Type") %in% colnames(probeTypes))) {
    stop("'probeTypes' must be a data.frame with a column 'Name' of probe ",
         "IDs and a column 'Type' indicating their design type")
  }
  
  # Construct 1-column matrix of beta values
  if (is(data, "MethylSet")) {
    betas <- as.matrix(getBeta(data))
  } else if (is.vector(data)) {
    betas <- matrix(
      data = data,
      ncol = 1,
      dimnames = list(names(data), NULL))
  } else {
    betas <- as.matrix(data)
  }
  
  # Compute densities
  betas <- betas[!is.na(betas), , drop = FALSE]
  n_probes <- nrow(betas)
  n_type1_grn <- sum(probeTypes$Type == "I-grn" &
                   probeTypes$Name %in% rownames(betas))
  n_type1_red <- sum(probeTypes$Type == "I-red" &
                       probeTypes$Name %in% rownames(betas))
  n_type2 <- sum(probeTypes$Type == "II" &
                   probeTypes$Name %in% rownames(betas))
  betas_density <- density(betas)
  
  typeI_grn_density <- suppressWarnings(
    density(
      x = betas[rownames(betas) %in%
                  probeTypes$Name[probeTypes$Type == "I-grn"]],
      weights = rep(1 / n_probes, n_type1_grn)))
  
  typeI_red_density <- suppressWarnings(
    density(
      x = betas[rownames(betas) %in%
                  probeTypes$Name[probeTypes$Type == "I-red"]],
      weights = rep(1 / n_probes, n_type1_red)))
  
  typeII_density <- suppressWarnings(
    density(
      x = betas[rownames(betas) %in%
                  probeTypes$Name[probeTypes$Type == "II"]],
      weights = rep(1 / n_probes, n_type2)))
  
  # Plot densities
  plot(
    x = betas_density,
    main = main,
    xlab = "Beta values",
    ylim = c(0, max(betas_density$y)),
    lwd = lwd,
    col = colors[1])
  lines(
    x = typeI_grn_density,
    col = colors[3],
    lty = 2,
    lwd = lwd)
  lines(
    x = typeI_red_density,
    col = colors[2],
    lty = 2,
    lwd = lwd)
  lines(
    x = typeII_density,
    col = colors[4],
    lty = 2,
    lwd = lwd)
  legend(
    x = legendPos,
    y = c("All probes", "Infinium I-grn", "Infinium I-red", "Infinium II"),
    col = colors,
    lwd = lwd,
    bg = "white",
    cex = cex.legend,
    lty = c(1, 2, 2))
}


## Plot betas with red/grn probes ------------------------------------------

par(mfrow = c(2,3))
for (i in rand_samp) {
  plotBetasByTypeColor(mset_noob[,i], 
                  main = paste0(i, " Noob"))
}
for (j in rand_samp) {
  plotBetasByTypeColor(mset_noob_swan[,j],
                  main = paste0(j, " Noob + SWAN"))
}


# Simulate qqplot from SeSAMe ---------------------------------------------

plotRedGrnQQ <- function(mset, main="R-G QQ Plot") {
  # Get Green Probes
  typeI_grn_M <- subsetProbes(mset, type = 'I-green', allele = "M")
  typeI_grn_U <- subsetProbes(mset, type = 'I-green', allele = "U")
  #Get Red Probes
  typeI_red_M <- subsetProbes(mset, type = 'I-red', allele = "M")
  typeI_red_U <- subsetProbes(mset, type = 'I-red', allele = "U")
  
  m <- max(c(typeI_grn_M, typeI_grn_U, typeI_red_M, typeI_red_U), na.rm=TRUE)
  
  qqplot(c(typeI_grn_M, typeI_grn_U), c(typeI_red_M, typeI_red_U),
    xlab = 'Infinium-I Red Signal', ylab = 'Infinium-I Grn Signal',
    main = main, cex = 0.5,
    xlim = c(0,m), ylim = c(0,m))
  graphics::abline(0,1,lty = 'dashed')
}

plot_redgrnqq <- function(raw_dat, noob_dat, swan_dat){
  
  par(mfrow=c(3,3), mar=c(3,3,2,1))
  
  #Make Raw Plots
  for (h in rand_samp){
    plotRedGrnQQ(raw_dat[,h], main= paste0(h, " (Raw)"))
    
    
  }
  #Make noob Plots
  for (i in rand_samp){
    plotRedGrnQQ(noob_dat[,i],
                          main= paste0(i, " (noob)"))
  }
  
  #Make noob + BMIQ Plots
  for (j in rand_samp){
    plotRedGrnQQ(swan_dat[,j],
                          main= paste0(j, " (noob + SWAN)"))
  }
}

plot_redgrnqq(raw_dat = mset_raw,
              noob_dat = mset_noob,
              swan_dat = mset_noob_swan)
# MDS Plots ---------------------------------------------------------------

#Allow MDS to plot more MDS coordinates
mdsPlot_n <- function(dat, numPositions = 1000, sampNames = NULL,
                      sampGroups = NULL, xlim, ylim, pch = 1,
                      pal = brewer.pal(8, "Dark2"), legendPos = "bottomleft",
                      legendNCol, main = NULL, n_mds = 2) {
  # Check inputs
  if (is(dat, "MethylSet") || is(dat, "RGChannelSet")) {
    b <- getBeta(dat)
  } else if (is(dat, "matrix")) {
    b <- dat
  } else {
    stop("dat must be an 'MethylSet', 'RGChannelSet', or 'matrix'.")
  }
  if (is.null(main)) {
    main <- sprintf(
      "Beta MDS\n%d most variable positions",
      numPositions)
  }
  
  o <- order(rowVars(b), decreasing = TRUE)[seq_len(numPositions)]
  d <- dist(t(b[o, ]))
  fit <- cmdscale(d, k = n_mds)
  if (missing(xlim)) xlim <- range(fit[, n_mds - 1]) * 1.2
  if (missing(ylim)) ylim <- range(fit[, n_mds]) * 1.2
  if (is.null(sampGroups)) sampGroups <- rep(1, numPositions)
  sampGroups <- as.factor(sampGroups)
  col <- pal[sampGroups]
  if (is.null(sampNames)) {
    plot(
      x = fit[, n_mds - 1],
      y = fit[, n_mds],
      col = col,
      pch = pch,
      xlim = xlim,
      ylim = ylim,
      xlab = paste0("MDS ", as.character(n_mds - 1)),
      ylab = paste0("MDS ", as.character(n_mds)),
      main = main)
  } else {
    plot(
      x = 0,
      y = 0,
      type = "n",
      xlim = xlim,
      ylim = ylim,
      xlab = paste0("MDS ", as.character(n_mds - 1)),
      ylab = paste0("MDS ", as.character(n_mds)),
      main = main)
    text(x = fit[, n_mds - 1], y = fit[, n_mds], sampNames, col = col)
  }
  numGroups <- length(levels(sampGroups))
  if (missing(legendNCol)) legendNCol <- numGroups
  if (numGroups > 1) {
    legend(
      x = legendPos,
      legend = levels(sampGroups),
      ncol = legendNCol,
      text.col = pal[seq_len(numGroups)])
  }
}

plotCpg <- function(dat, cpg, pheno, type = c("categorical", "continuous"),
                    measure = c("beta", "M"), ylim = NULL, ylab = NULL,
                    xlab = "", fitLine = TRUE, mainPrefix = NULL,
                    mainSuffix = NULL) {
  if (is.numeric(cpg)) cpg <- rownames(dat)[cpg]
  type <- match.arg(type)
  measure <- match.arg(measure)
  if (is(dat, "MethylSet") || is(dat, "RGChannelSet")) {
    if (measure == "beta") {
      # NOTE: as.matrix() necessary in case 'dat' is a
      #       DelayedArray-backed minfi object
      x <- as.matrix(getBeta(dat[cpg, ]))
      if (is.null(ylab)) ylab <- "Beta"
      if (is.null(ylim)) ylim <- c(0, 1)
    } else if (measure == "M") {
      x <- getM(dat[cpg, ])
      if (is.null(ylab)) ylab <- "M"
      if (is.null(ylim)) ylim <- range(x)
    }
  } else {
    if (is.null(ylab)) ylab <- "unknown"
    x <- dat
    if (is.vector(x)) {
      x <- matrix(x, ncol = 1)
    } else {
      x <- as.matrix(x)
    }
  }
  main <- paste(mainPrefix, cpg, mainSuffix)
  names(main) <- cpg
  for (i in cpg) {
    if (type == "categorical") {
      stripchart(
        x = x[i, ] ~ pheno,
        vertical = TRUE,
        method = "jitter",
        jitter = 0.15,
        ylab = ylab,
        ylim = ylim,
        main = main[i])
    } else if (type == "continuous") {
      plot(
        x = pheno,
        y = x[i,],
        ylab = ylab,
        xlab = xlab,
        ylim = ylim,
        main = main[i])
      if (fitLine) abline(lm(x[i, ] ~ pheno), col = "blue")
    }
  }
}
#Code to prep data and make mds plots
make_mds_plots_n <- function(vals, meta, plot = c("sex", "vape", "center"), nump = 1000, lim = c(-10, 10), mlim = c(-50, 50), n_mds = 2){
  #Get only samples with methylation data and arrange in same order as betas so that you can use for mds plots
  meta <- meta %>% 
    filter(sentrix_name %in% colnames(vals)) %>% 
    arrange(factor(sentrix_name, levels = colnames(vals)))
  
  #Get mvals
  mvals <- BetaValueToMValue(vals)
  
  # Set Color Pallattes
  sex_pal <- c("deepskyblue", "deeppink3")
  center_colors <- c("cyan3", "chartreuse3", "darkorange")
  
  
  #MDSPlots (Betas)
  par(mfrow = c(1,2))
  if ("sex" %in% plot) {
    #Clustering by Sex
    mdsPlot_n(vals, sampGroups = meta$sex_lab, sampNames = meta$sid,
              numPositions=nump, ylim = lim, main = paste0("Sex - Betas - ", as.character(nump)), 
              pal = sex_pal,
              n_mds = n_mds)
    
    mdsPlot_n(mvals, sampGroups = meta$sex_lab, sampNames = meta$sid,
              numPositions=nump, ylim = mlim, main = paste0("Sex - M-values - ", as.character(nump)), 
              pal = sex_pal,
              n_mds = n_mds)
  }
  if ("vape" %in% plot) {
    #Clustering by Vape Status
    mdsPlot_n(vals, sampGroups = meta$vape_text, sampNames = meta$sid,
              numPositions=nump, ylim = lim, main = paste0("Vape Status - Betas - ", as.character(nump)),
              n_mds = n_mds)
    
    mdsPlot_n(mvals, sampGroups = meta$vape_text, sampNames = meta$sid,
              numPositions=nump, ylim = mlim, main = paste0("Vape Status - M-values - ", as.character(nump)),
              n_mds = n_mds)
  }
  if ("center" %in% plot) {
    #Clustering by center
    mdsPlot_n(vals, sampGroups = meta$recruitment_center, sampNames = meta$sid,
              numPositions=nump, 
              ylim = lim, 
              pal = center_colors, 
              main = paste0("Recruitment Center - Betas - ", as.character(nump)),
              n_mds = n_mds)
    
    mdsPlot_n(mvals, sampGroups = meta$recruitment_center, sampNames = meta$sid,
              numPositions=nump, ylim = mlim, pal = center_colors, 
              main = paste0("Recruitment Center - M-values - ", as.character(nump)),
              n_mds = n_mds)
  }
}

make_mds_plots_n(ns_betas, metadata_sex, plot = "vape", nump = 1000, n_mds = 2)
make_mds_plots_n(ns_betas, metadata_sex, plot = "vape", nump = 10000, lim = c(-10, 10), mlim = c(-100, 100), n_mds = 2)
