library(sesameData)
library(tidyverse)
library(data.table)

methyl_betas <- fread("./Output/2022_07_27_MethylationDat/betas.tsv")
methyl_mvals <- fread("./Output/2022_07_27_MethylationDat/mvals.tsv")

# sesame betas sorted by probenames
# the manifest sorted by chr -> position
EPIC.hg38.manifest <- sesameDataGet('EPIC.hg38.manifest')
EPIC.hg38.manifest_df <-
  EPIC.hg38.manifest %>%
  as_tibble() %>%
  bind_cols(probe = EPIC.hg38.manifest@ranges@NAMES, .) %>%
  arrange(probe)

identical(methyl_betas$Probe, EPIC.hg38.manifest_df$probe)

filter_probes_autosomal <-
  EPIC.hg38.manifest_df$seqnames %in% paste0("chr", 1:22)
filter_probes_nonSNP <- 
  !EPIC.hg38.manifest_df$probeType == "rs"


filter_probes_anyNA <-
  methyl_betas %>% apply(., 1, function(x) { is.na(x) %>% any})

filter_probes_autosomal <-
  EPIC.hg38.manifest_df$seqnames %in% paste0("chr", 1:22)

probes_beta_range <-
  apply(methyl_betas[ , -1], 1, function(x) { abs(max(x) - min(x)) })
filter_beta_range <- 
  # probes_beta_range > 0.05
  order(probes_beta_range) %in% c(1:10000)

filter_probes <- !filter_probes_anyNA & filter_beta_range & filter_probes_autosomal






pca_results <-
  prcomp(methyl_mvals[filter_probes, -1] %>% t, center = T, scale = F)
pca_percent_exp <-
  (pca_results$sdev^2) %>% `/`(., sum(.)) %>% `*`(100) %>% formatC(digits = 2) %>%
  gsub(" ", "", .) %>%
  paste0(" (", ., "%)") %>%
  paste0("PC", 1:length(.), .)



pca_df <-
  metadata_joined %>%
  bind_cols(pca_results$x[ , 1:3] %>% as_tibble)

ggplot(pca_df, aes(PC1, PC2, color = sex_lab, shape = sex_lab)) + #, color = sex, shape = sex)) +
  geom_point(alpha = 0.75) +
  geom_text_repel(
  aes(label = new_id),
  color = "black", size = 2) +
  # scale_shape_manual(name = "Sex", values = palette_shape_sex, labels = c("Female", "Male")) +
  # scale_color_manual(name = "Sex", values = palette_color_sex, labels = c("Female", "Male")) +
  xlab(pca_percent_exp[1]) + ylab(pca_percent_exp[2]) +
  ggtitle("Methylation PCA") +
  theme_few() +
  facet_grid(~ recruitment_center) +
  geom_blank()





