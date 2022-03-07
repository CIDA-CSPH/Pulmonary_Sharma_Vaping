# 60651 x 49
# convert the IDs given to the core to PIDs

run_DESeq_lung_fxn <-
  function(variable_of_interest, contrast_extract = NULL,
         suffix = "_RUV",
         outlier_downweight = F) {
  
metadata_IOS <-
  metadata_naive %>%
  filter(SID != "102") %>%
  mutate(Feature := !!sym(variable_of_interest)) %>%
  filter(!is.na(Feature)) %>%
  mutate(AgeScaled = scale(Age),
         GenderMale = Gender == "1")

vstcounts_preruv <-
  DESeqDataSetFromMatrix(
    countData =
      counts_mat_naive %>%
      .[ , colnames(.) %in% metadata_IOS$NewID],
    colData = metadata_IOS,
    design = ~ VapingLast6Mo + AgeScaled + GenderMale + EthnicityLatino
  ) %>%
  DESeq() %>%
  vst %>%
  assay

resid_vals <-
  vstcounts_preruv %>%
  apply(., 1, function(y) {
    lm(y ~ Feature +
         AgeScaled + EthnicityLatino + GenderMale,
       metadata_IOS) %>% resid
  })

ruv_results <-
  RUVr(x = vstcounts_preruv,
       k = 4,
       residuals = resid_vals %>% t,
       isLog = T)

metadata_IOS <-
  metadata_IOS %>%
  bind_cols(ruv_results$W %>% as_tibble() %>% set_names(paste0("RUV", 1:ncol(.))))

PCA_ruv <-
  prcomp(t(ruv_results$normalizedCounts))$x
ggplot(metadata_IOS %>%
         bind_cols(PCA_ruv[ , 1:2] %>% as_tibble),
       aes(PC1, PC2, color = RecruitmentCenter, shape = VapingLast6Mo)) +
  geom_point() +
  scale_shape_manual(values = palette_shape_vaping) +
  theme_few()

deseqobj <-
      DESeqDataSetFromMatrix(
        countData = counts_mat_naive %>%
          .[ , colnames(.) %in% metadata_IOS$NewID],
        colData = metadata_IOS,
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
      assays(deseqobj)[["weights"]][cooks > 5] <- 1e-5
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
    
    
    assign(x = paste0("ruv_results_", variable_of_interest, suffix),
           value = ruv_results, envir = .GlobalEnv)
    assign(x = paste0("deseqobj_", variable_of_interest, suffix),
           value = deseqobj, envir = .GlobalEnv)
    assign(x = paste0("DEResults_", variable_of_interest, suffix),
           value = DEResults_out, envir = .GlobalEnv)
    
  }


list_genes_DE_vaping <-
  DEResults_VapingLast6Mo_RUV %>%
  filter(padj < 0.05) %>%
  .$row

sapply(c("AX", "Fres", "R5", "R20",
         "R20R5", "X5", "X20", "X20X5", "VT"),
       run_DESeq_lung_fxn, outlier_downweight = T)


# results(deseqobj_X20_RUV, name = "Feature", tidy = T) %>%
#   arrange(padj) %>%
#   left_join(., gene_annotations_short, by = c("row" = "ENSG")) %>% head
