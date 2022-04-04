# run_DESeq_lung_fxn <-
#   function(variable_of_interest,
#            contrast_extract = NULL,
#            outputname = paste0("DEResults_", variable_of_interest, suffix),
#            suffix = "_ageethn_noaurora",
#            outlier_downweight = F) {
#     
#     metadata_final <-
#       metadata_final %>%
#       mutate(Feature := !!sym(variable_of_interest))
# 
#     filter_samples <-
#       (metadata_final$SID != "102") &
#         # (metadata_final$RecruitmentCenter == "Pueblo") &
#         (metadata_final$RecruitmentCenter != "Aurora") &
#         (!is.na(metadata_final$Feature))
# 
#     final_counts_mat <-
#       count_mat[filter_genes, c(T, filter_samples)]
# 
#     final_counts_mat <-
#       final_counts_mat[, -1] %>%
#       as.matrix() %>%
#       `rownames<-`(., final_counts_mat$Feature)
# 
#     metadata_final <-
#       metadata_final %>%
#       filter(filter_samples)
# 
#     deseqobj <-
#       DESeqDataSetFromMatrix(
#         countData = final_counts_mat,
#         colData = metadata_final,
#         # design = ~ Feature
#         # design = ~ Feature + Age
#         design = ~ Feature + Age + Latino + GenderMale
#         # design = ~ Feature + RecruitmentCenter
#         # design = ~ Feature + CityAurora
#       )
# 
#     deseqobj <-
#       deseqobj %>%
#       DESeq()
# 
#     if (outlier_downweight) {
#       # https://support.bioconductor.org/p/103959/
#       cooks <- assays(deseqobj)[["cooks"]]
#       assays(deseqobj, withDimnames = F)[["weights"]] <-
#         matrix(1, nrow(deseqobj), ncol(deseqobj))
#       assays(deseqobj)[["weights"]][cooks > 10] <- 1e-5
#       deseqobj <-
#         deseqobj %>%
#         DESeq()
#     }
# 
#     if (!is.null(contrast_extract)) {
#       DEResults_out <-
#         results(deseqobj, contrast = contrast_extract, tidy = T) %>%
#         arrange(padj) %>%
#         left_join(., gene_annotations_short, by = c("row" = "ENSG"))
#     } else {
#       DEResults_out <-
#         results(deseqobj, name = "Feature", tidy = T) %>%
#         arrange(padj) %>%
#         left_join(., gene_annotations_short, by = c("row" = "ENSG"))
#     }
# 
#     assign(x = outputname, value = DEResults_out, envir = .GlobalEnv)
#   }
# 
# 
# 
# 
# # -- replacing outliers and refitting for 693 genes
# metadata_naive$VapingLast6Mo %>% table()
# run_DESeq_lung_fxn("VapingLast6Mo",
#   contrast_extract = c("Feature", "Vaping", "Control")
# )
# DEResults_VapingLast6Mo %>%
#   filter(padj < 0.1) %>%
#   head()
# run_DESeq_lung_fxn("VapingTimesPer",
#   contrast_extract = c("Feature", "Vaping", "Control")
# )
# DEResults_VapingTimesPer %>%
#   filter(padj < 0.1) %>%
#   head()
# run_DESeq_lung_fxn("VapingPast30",
#   contrast_extract = c("Feature", "Vaping", "Control")
# )
# DEResults_VapingPast30 %>%
#   filter(padj < 0.05) %>%
#   head()
# 
# 
# run_DESeq_lung_fxn("AX", outlier_downweight = T)
# DEResults_AX %>% filter(padj < 0.1)
# run_DESeq_lung_fxn("VT", outlier_downweight = T)
# DEResults_VT %>% filter(padj < 0.1)
# run_DESeq_lung_fxn("Fres", outlier_downweight = T)
# DEResults_Fres %>% filter(padj < 0.1)
# run_DESeq_lung_fxn("R5", outlier_downweight = T)
# DEResults_R5 %>% filter(padj < 0.1)
# run_DESeq_lung_fxn("X5", outlier_downweight = T)
# DEResults_X5 %>% filter(padj < 0.1)
# 
# 
# 
# 
# 
# plotCounts(deseqobj_naive[, metadata_final$RecruitmentCenter == "Pueblo"],
#   gene = "ENSG00000000938.13", intgroup = "VapingLast6Mo"
# )
# 
# 
# 
# 
# DEResults_VapingLast6Mo_pueblocenter %>% .[1:10, ]
# plotCounts(deseqobj_naive[, metadata_final$RecruitmentCenter == "Pueblo"],
#   gene = DEResults_VapingLast6Mo_pueblocenter$row[2], intgroup = "VapingLast6Mo"
# )
# DEResults_VapingLast6Mo_ageethn_noaurora %>% .[1:10, ]
# plotCounts(deseqobj_naive[, metadata_final$RecruitmentCenter != "Aurora"],
#   gene = DEResults_VapingLast6Mo_ageethn_noaurora$row[1], intgroup = "VapingLast6Mo"
# )










deseqobj_VapingLast6Mo_naive <-
  DESeqDataSetFromMatrix(
    countData =
      counts_mat_naive %>%
      .[ , colnames(.) %in% metadata_final$NewID],
    colData = metadata_final,
    design = ~ Feature + Age + Latino + GenderMale
  ) %>%
  DESeq()
