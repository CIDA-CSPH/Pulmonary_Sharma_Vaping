# 60651 x 49
# convert the IDs given to the core to PIDs
metadata_final <-
  metadata_naive %>%
  filter(SID != "102") %>%
  mutate(AgeScaled = scale(Age))
 
vstcounts_preruv <-
  DESeqDataSetFromMatrix(
    countData =
      counts_mat_naive %>%
      .[ , colnames(.) %in% metadata_final$NewID],
    colData = metadata_final,
    design = ~ VapingLast6Mo + AgeScaled + GenderMale + EthnicityLatino
) %>%
  DESeq() %>%
  vst %>%
  assay

resid_vals <-
  vstcounts_preruv %>%
  apply(., 1, function(y) {
    lm(y ~ VapingLast6Mo +
         AgeScaled + EthnicityLatino + GenderMale,
       metadata_final) %>% resid
  })

dm_ruv_test <-
  vstcounts_preruv %>%
  .[filter_genes_variance, ] %>%
  t %>% dist

# ruv_out <-
#   sapply(1:10,
#          function(ktry) {
#            ruv_results <-
#              RUVr(x = vstcounts_preruv,
#                   k = ktry,
#                   residuals = resid_vals %>% t,
#                   isLog = T)
#            
#            adonis(dm_ruv_test ~
#                     ruv_results$W +
#                     VapingLast6Mo + AgeScaled + EthnicityLatino + GenderMale,
#                   data = metadata_final) %>% .$aov.tab %>% .$R2 %>% .[6]
#          }
#   )
# 
# ggplot(data = NULL, aes(1:10, 1 - ruv_out)) +
#   geom_point() +
#   xlab("# of RUVr Components") +
#   ylab("% Expression Variance Explained\n(by all covariates)") +
#   theme_few()
# 


ruv_results <-
  RUVr(x = vstcounts_preruv,
       k = 4,
       residuals = resid_vals %>% t,
       isLog = T)
adonis(dm_ruv_test ~
         ruv_results$W +
         VapingLast6Mo + AgeScaled + EthnicityLatino + GenderMale,
       data = metadata_final)

metadata_final <-
  metadata_final %>%
  bind_cols(ruv_results$W %>% as_tibble() %>% set_names(paste0("RUV", 1:ncol(.))))



run_DESeq_lung_fxn <-
  function(variable_of_interest, contrast_extract = NULL,
           suffix = "_RUV",
           outlier_downweight = F) {
    
    metadata_final <-
      metadata_final %>%
      mutate(Feature := !!sym(variable_of_interest))
  
    deseqobj <-
      DESeqDataSetFromMatrix(
        countData = counts_mat_naive %>%
          .[ , colnames(.) %in% metadata_final$NewID],
        colData = metadata_final,
        design = ~ Feature +
          RUV1 + RUV2 + RUV3 + RUV4 +
          AgeScaled + EthnicityLatino + GenderMale
      )

    deseqobj <-
      deseqobj %>%
      DESeq()
    
    if (outlier_downweight) {
      # https://support.bioconductor.org/p/103959/
      cooks <- assays(deseqobj)[["cooks"]]
      assays(deseqobj, withDimnames = F)[["weights"]] <-
        matrix(1, nrow(deseqobj), ncol(deseqobj))
      assays(deseqobj)[["weights"]][cooks > 10] <- 1e-5
      deseqobj <-
        deseqobj %>%
        DESeq()
    }
    
    if (!is.null(contrast_extract)) {
      DEResults_out <-
        results(deseqobj, contrast = contrast_extract, tidy = T) %>%
        arrange(padj) %>%
        left_join(., gene_annotations_short, by = c("row" = "ENSG"))
    } else {
      DEResults_out <-
        results(deseqobj, name = "Feature", tidy = T) %>%
        arrange(padj) %>%
        left_join(., gene_annotations_short, by = c("row" = "ENSG"))
    }

    assign(x = paste0("deseqobj_", variable_of_interest, suffix),
           value = deseqobj, envir = .GlobalEnv)
    assign(x = paste0("DEResults_", variable_of_interest, suffix),
           value = DEResults_out, envir = .GlobalEnv)
  }



# -- replacing outliers and refitting for 693 genes
metadata_naive$VapingLast6Mo %>% table
run_DESeq_lung_fxn("VapingLast6Mo",
                   contrast_extract = c("Feature", "Vaping", "Control")
)

DEResults_VapingLast6Mo_RUV %>%
  filter(padj < 0.05)

plotCounts(deseqobj_VapingLast6Mo_RUV,
           gene = DEResults_VapingLast6Mo_RUV$row[4],
           intgroup = "VapingLast6Mo",
           returnData = T) %>%
  ggplot(data = ., aes(VapingLast6Mo, count)) +
  geom_quasirandom(width = 0.10)

ggplot(data = NULL, aes(metadata_final$VapingLast6Mo, 
                     ruv_results$normalizedCounts[
                       DEResults_VapingLast6Mo_RUV$row[4], ])) +
  geom_quasirandom(width = 0.10) 






PCA_ruv <-
  prcomp(t(ruv_results$normalizedCounts), scale. = T)
percentexp_ruv <-
  formatC(PCA_ruv$sdev^2 / sum(PCA_ruv$sdev^2) * 100, digits = 1, format = "f") %>%
  paste0("PC", 1:length(.), " (", ., "%)")

percentexp_ruv <-
  formatC(PCA_ruv$sdev^2 / sum(PCA_ruv$sdev^2) * 100, digits = 1, format = "f") %>%
  paste0("PC", 1:length(.), " (", ., "%)")

metadata_final %>% 
  bind_cols(PC1 = PCA_ruv$x[ , 1],
            PC2 = PCA_ruv$x[ , 2]) %>%
  ggplot(.,
       aes(PC1, PC2,
           shape = VapingLast6Mo,
           color = RecruitmentCenter)) +
  geom_point(size = 4, alpha = 0.75) +
  theme_few() + 
  xlab(percentexp_ruv[1]) + ylab(percentexp_ruv[2]) +
  scale_color_manual("Center", values = palette_color_center) +
  scale_shape_manual("Exposure", values = palette_shape_vaping) +
  NULL




metadata_naive %>% 
  bind_cols(PC1 = pca_naive$x[ , 1],
            PC2 = pca_naive$x[ , 2]) %>%
  filter(SID != "102") %>%
  ggplot(.,
         aes(PC1, PC2,
             shape = VapingLast6Mo,
             color = RecruitmentCenter)) +
  geom_point(size = 4, alpha = 0.75) +
  theme_few() + 
  xlab(percentexp_naive[1]) + ylab(percentexp_naive[2]) +
  scale_color_manual("Center", values = palette_color_center) +
  scale_shape_manual("Exposure", values = palette_shape_vaping) +
  NULL


plot_grid(
  ggplot(metadata_all,
       aes(RecruitmentCenter, `260/280`,
           color = RecruitmentCenter)) +
  geom_quasirandom(width = 0.2) +
  theme_few() +
  scale_color_manual(values = palette_color_center) +
  stat_summary(color = "black") +
  stat_summary(color = "black", geom = "errorbar", width = 0.1) +
  theme(legend.position = "none") +
    xlab("Center"),
ggplot(metadata_all,
       aes(RecruitmentCenter, `Total RNA (ng)`,
           color = RecruitmentCenter)) +
  geom_quasirandom(width = 0.2) +
  theme_few() +
  stat_summary(color = "black") +
  stat_summary(color = "black", geom = "errorbar", width = 0.1) +
  scale_color_manual(values = palette_color_center) +
  theme(legend.position = "none") +
  xlab("Center"),
nrow = 2
)

boxplot(metadata_all$`Total RNA (ng)` ~
          metadata_all$RecruitmentCenter)











