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
M <- mvals %>% select(-CpG_Site, -targets$sentrix_name[is.na(targets$vape_6mo_lab)])

### Make sure clinical metadata same arrangement as Mvals
clin_metadata <- clin_metadata[match(names(M), clin_metadata$sentrix_name),]
# Drop the NA Vaper From Targets-------------------------------------------------------
targets <- targets %>% drop_na(vape_6mo_lab)

# setup the factor of interest
grp <- factor(clin_metadata$vape_6mo_lab, levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"), labels = c(0,1))

z <- data.frame(#sex_lab = factor(clin_metadata$sex_lab, levels = c("Female", "Male"), labels = c(0,1)), 
                age = clin_metadata$age)
rgSet <- read.metharray.exp(targets = targets)

# extract Illumina negative control data
INCs <- getINCs(rgSet)
head(INCs)

# add negative control data to M-values
Mc <- rbind(M,INCs)

# create vector marking negative controls in data matrix
ctl1 <- rownames(Mc) %in% rownames(INCs)

rfit1 <- RUVfit(Y = Mc, X = grp, Z = z, ctl = ctl1) # Stage 1 analysis

rfit2 <- RUVadj(Y = Mc, fit = rfit1)

top1 <- topRUV(rfit2, num=Inf, p.BH = 1)
head(top1)

ctl2 <- rownames(M) %in% rownames(top1[top1$p.BH_X1.1 > 0.5,])
table(ctl1)

# Perform RUV adjustment and fit
rfit3 <- RUVfit(Y = M, X = grp, Z = z, ctl = ctl2) # Stage 2 analysis
rfit4 <- RUVadj(Y = M, fit = rfit3)

# Look at table of top results
topRUV(rfit4)
table(ctl2)

#Get adjusted values
Madj <- getAdj(M, rfit3)

vape_cols <- if_else(grp == 0, brewer.pal(3, "Dark2")[1], brewer.pal(3, "Dark2")[2])

par(mfrow=c(1,2))
plotMDS(M, labels=clin_metadata$sid, col= vape_cols,
        main="Unadjusted", gene.selection = "common")
legend("topright",legend=c("No Vape", "Vape"),pch=16,cex=1,col=unique(vape_cols))
plotMDS(Madj, labels=clin_metadata$sid, col= vape_cols,
        main="Adjusted: RUV-inverse", gene.selection = "common", xlim = c(-1, 1), ylim = c(-1, 1))
legend("topright",legend=c("No Vape", "Vape"),pch=16,cex=1,col=unique(vape_cols))
