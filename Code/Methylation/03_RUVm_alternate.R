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
library(EDASeq)

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

# setup the factors of interest
vape_status <- factor(clin_metadata$vape_6mo_lab, levels = c("Did Not Vape in Last 6 Months", "Vaped in Last 6 Months"), labels = c(0,1))

sex <- factor(clin_metadata$sex_lab, levels = c("Female", "Male"), labels = c(0,1))

age <- scale(clin_metadata$age, scale = T, center = F)

center <- factor(clin_metadata$recruitment_center, levels = c("Aurora", "Pueblo", "CommCity/Denver")) #for plotting only

#Design Matrix
des <- model.matrix(~vape_status + sex + age)
des

# limma differential methylation analysis
lfit1 <- lmFit(M, design=des)
lfit2 <- eBayes(lfit1) # Stage 1 analysis
# Look at table of top results
topTable(lfit2)

#use to define ECPs
topl1 <- topTable(lfit2, num=Inf)

head(topl1)

ctl3 <- rownames(M) %in% rownames(topl1[topl1$adj.P.Val > 0.5,])
table(ctl3)

# Perform RUV adjustment and fit
rfit5 <- RUVfit(Y = M, X = vape_status, Z = data.frame("sex" = sex, "age" = age), ctl = ctl3) # Stage 2 analysis
rfit6 <- RUVadj(Y = M, fit = rfit5)
# Look at table of top results
topRUV(rfit6)
#Get adjusted values
Madj <- getAdj(M, rfit5)

# Use RUV-4 in stage 2 of RUVm with k=1 and k=2
rfit7 <- RUVfit(Y = M, X = vape_status, Z = data.frame("sex" = sex, "age" = age), ctl = ctl3,
                method = "ruv4", k=1) # Stage 2 with RUV-4, k=1
rfit9 <- RUVfit(Y = M, X = vape_status, Z = data.frame("sex" = sex, "age" = age), ctl = ctl3,
                method = "ruv4", k=2) # Stage 2 with RUV-4, k=1
# get adjusted values
Madj1 <- getAdj(M, rfit7)
Madj2 <- getAdj(M, rfit9)

#set colors
vape_cols <- if_else(vape_status == 0, brewer.pal(3, "Dark2")[1], brewer.pal(3, "Dark2")[2])
recruit_cols <- case_when(center == "Aurora" ~ brewer.pal(3, "Set1")[1],
                          center == "Pueblo" ~ brewer.pal(3, "Set1")[2],
                          center == "CommCity/Denver" ~ brewer.pal(3, "Set1")[3])

# Vape Plots
par(mfrow=c(2,2))

plotMDS(M, labels=targets$Sample_Name, col=vape_cols,
        main="Unadjusted", gene.selection = "common")
legend("top",legend=c("Not Vaped","Vaped"),pch=16,cex=1,col=unique(vape_cols))

plotMDS(Madj, labels=targets$Sample_Name, col=vape_cols,
        main="Adjusted: RUV-inverse", gene.selection = "common")
legend("topright",legend=c("Not Vaped","Vaped"),pch=16,cex=1,col=unique(vape_cols))

plotMDS(Madj1, labels=targets$Sample_Name, col=vape_cols,
        main="Adjusted: RUV-4, k=1", gene.selection = "common")
legend("bottom",legend=c("Not Vaped","Vaped"),pch=16,cex=1,col=unique(vape_cols))

plotMDS(Madj2, labels=targets$Sample_Name, col=vape_cols,
        main="Adjusted: RUV-4, k=2", gene.selection = "common")
legend("bottomright",legend=c("Not Vaped","Vaped"),pch=16,cex=1,col=unique(vape_cols))

#Recruitment Center Plots
par(mfrow=c(2,2))

plotMDS(M, labels=targets$Sample_Name, col=recruit_cols,
        main="Unadjusted", gene.selection = "common")
legend("top",legend=c("Aurora","Pueblo","CommCity/Denver"),pch=16,cex=1,col=unique(recruit_cols))

plotMDS(Madj, labels=targets$Sample_Name, col=recruit_cols,
        main="Adjusted: RUV-inverse", gene.selection = "common")
legend("topright",legend=c("Aurora","Pueblo","CommCity/Denver"),pch=16,cex=1,col=unique(recruit_cols))

plotMDS(Madj1, labels=targets$Sample_Name, col=recruit_cols,
        main="Adjusted: RUV-4, k=1", gene.selection = "common")
legend("bottom",legend=c("Aurora","Pueblo","CommCity/Denver"),pch=16,cex=1,col=unique(recruit_cols))

plotMDS(Madj2, labels=targets$Sample_Name, col=recruit_cols,
        main="Adjusted: RUV-4, k=2", gene.selection = "common")
legend("bottomright",legend=c("Aurora","Pueblo","CommCity/Denver"),pch=16,cex=1,col=unique(recruit_cols))


# RLE Plots ---------------------------------------------------------------

rle_k0 <- ruv::ruv_rle(M %>% t() %>% as.matrix(), ylim = c(-1,1), 
                       rowinfo = clin_metadata$vape_6mo_lab) + 
  theme(axis.title.y = element_text(angle = 90, size = 12)) + 
  labs(y = "K = 0") + 
  geom_hline(yintercept = c(-0.25, 0.25), col = 'red', linetype = "dashed") +
  scale_fill_manual(values = unique(vape_cols))

rle_k1 <- ruv::ruv_rle(Madj1 %>% t() %>% as.matrix(), ylim = c(-1,1), 
                       rowinfo = clin_metadata$vape_6mo_lab) + 
  theme(axis.title.y = element_text(angle = 90, size = 12)) + 
  labs(y = "K = 1") +
  geom_hline(yintercept = c(-0.25, 0.25), col = 'red', linetype = "dashed")+
  scale_fill_manual(values = unique(vape_cols))

rle_k2 <- ruv::ruv_rle(Madj2 %>% t() %>% as.matrix(), ylim = c(-1,1), 
                       rowinfo = clin_metadata$vape_6mo_lab) +
  theme(axis.title.y = element_text(angle = 90, size = 12)) + 
  labs(y = "K = 2", ) +
  geom_hline(yintercept = c(-0.25, 0.25), col = 'red', linetype = "dashed")+
  scale_fill_manual(values = unique(vape_cols))

ggpubr::ggarrange(rle_k0, rle_k1, rle_k2, nrow = 3, ncol = 1, 
                  legend = "bottom", common.legend = T)

