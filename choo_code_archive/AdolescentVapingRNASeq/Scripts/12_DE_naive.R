# replacing outliers and refitting for 2193 genes
metadata_naive$AgeScaled <-
  scale(metadata_naive$Age)
deseqobj_naive <-
  DESeqDataSetFromMatrix(countData = counts_mat_naive,
                         colData = metadata_naive,
                         design = ~ VapingLast6Mo + AgeScaled + EthnicityLatino + GenderMale
  ) %>%
  DESeq()

# plotDispEsts(deseqobj_naive)
# looks like she normalised using a VST transformation (Trent)
vstcounts_naive <-
  deseqobj_naive %>% vst()

#DERresults(deseqobj_naive, tidy = T)
results(deseqobj_naive, tidy = T) #Looks like old function might be depractaed -Trent

#more gene filtering? (Trent)

filter_genes_variance <- 
  vstcounts_naive %>%
  assay %>%
  apply(., 1, mean) %>%
  order(decreasing = T) %>%
  .[1:2500]

#PCA start (Trent)
pca_naive <-
  assay(vstcounts_naive) %>%
  .[filter_genes_variance, ] %>%
  t %>%
  prcomp(scale = T)
percentexp_naive <-
  formatC(pca_naive$sdev^2 / sum(pca_naive$sdev^2) * 100, digits = 1, format = "f") %>%
  paste0("PC", 1:length(.), " (", ., "%)")
metadata_for_pca_naive <-
  metadata_naive %>%  
  bind_cols(PC1 = pca_naive$x[ , 1],
            PC2 = pca_naive$x[ , 2])
ggplot(metadata_for_pca_naive,
       aes(PC1, PC2,
           # shape = VapingLast6Mo,
           color = RecruitmentCenter)) +
  geom_point(size = 4, alpha = 0.75) +
  geom_text_repel(aes(label = if_else(PC1 < -50, as.character(SID), "")),
                  show.legend = F) +
  theme_few() + theme(legend.position = "bottom") +
  xlab(percentexp_naive[1]) + ylab(percentexp_naive[2]) +
  scale_color_manual("Center", values = palette_color_center) +
  # scale_shape_manual("Group", values = palette_shape_vaping) +
  NULL

