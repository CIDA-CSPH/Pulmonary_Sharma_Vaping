

prep_DE_tables <-
  function(results_table,
           asterisk_for_vaping = F) {
    out <-
      results_table %>%
      transmute(ENSG = row,
                Symbol = symbol,
                log2FoldChange,
                `p-value` = pvalue,
                FDR = padj,
                Signif = if_else(FDR < 0.05, "*", ""))
    
    if (asterisk_for_vaping) {
      tmp <-
        tibble(row = results_table$row) %>%
        left_join(DEResults_VapingLast6Mo_RUV, by = "row") %>%
        mutate(direction = if_else(log2FoldChange < 0, "↓", "↑"),
               vaping_signif = if_else(padj < 0.05, direction, ""))
        
        
        out$`DE w/ Vaping` <- tmp$vaping_signif
    }
    
    out
  }


ls(pattern = "DEResults*") %>%
  sort


ls(pattern = "DEResults*") %>%
  sapply(., function(x) {
    get(x) %>% filter(padj < 0.05) %>% nrow
  }) %>%
  enframe() %>%
  transmute(Feature = gsub("DEResults_|_RUV", "", name),
            `# Signif (FDR < 0.05)` = value) %>%
  arrange(!grepl("Vaping", Feature)) %>%
  View()


ls(pattern = "DEResults*") %>%
  sapply(., function(x) {
    get(x) %>% filter(padj < 0.05 & row %in% list_genes_DE_vaping) %>% nrow
  }) %>%
  enframe() %>%
  transmute(Feature = gsub("DEResults_|_RUV", "", name),
            `# Also DE in Vaping` = value) %>%
  arrange(!grepl("Vaping", Feature)) %>%
  View()


write_xlsx(
  list(VapingLast6Mo = DEResults_VapingLast6Mo_RUV %>% prep_DE_tables(),
       GSEA_VapingLast6Mo = GSEAResults_Vaping,
       AX = DEResults_AX_RUV %>% prep_DE_tables(T),
       Fres = DEResults_Fres_RUV %>% prep_DE_tables(T),
       R5 = DEResults_R5_RUV %>% prep_DE_tables(T),
       R20 = DEResults_R20_RUV %>% prep_DE_tables(T),
       `R20-R5` = DEResults_R20R5_RUV %>% prep_DE_tables(T),
       VT = DEResults_VT_RUV %>% prep_DE_tables(T),
       X5 = DEResults_X5_RUV %>% prep_DE_tables(T),
       X20 = DEResults_X20_RUV %>% prep_DE_tables(T),
       `X20-X5` = DEResults_X20X5_RUV %>% prep_DE_tables(T)
       ) %>%
    set_names(., names(.) %>% paste0(1:length(.), "_", .)),
  path = "./Output/20210412_Vaping_RNAseq_v2.xlsx"
)
# 
# 
# 
# 
# library(eulerr)
# list(Last6Mo = DEResults_VapingLast6Mo_age_ethn,
#      Last30 = DEResults_VapingPast30_age_ethn,
#      TimesDay = DEResults_VapingTimesPer_age_ethn
# ) %>%
#   map(~ filter(.x, padj < 0.1) %>% .$row) %>%
#   euler() %>%
#   plot(quantities = T)
# 
# 
# set.seed(1234)
# list(Mod1 = DEResults_VapingLast6Mo,
#      Mod2 = DEResults_VapingLast6Mo_age,
#      Mod3 = DEResults_VapingLast6Mo_age_ethn,
#      Mod4 = DEResults_VapingLast6Mo_center
# ) %>%
#   map(~ filter(.x, padj < 0.1) %>% .$row) %>%
#   euler() %>%
#   plot(quantities = T)
# 
# 
# 
# index_feature <- 4
# tmp_exp <-
#   assay(vstcounts)[
#     DEResults_VapingLast6Mo$row[index_feature],
#     metadata_naive$SID != "102"]
# ggplot(data = metadata_naive %>% filter(SID != "102"),
#        aes(VapingLast6Mo, tmp_exp, color = VapingLast6Mo)) +
#   geom_quasirandom(alpha = 0.8, width = 0.2) +
#   theme_few() +
#   xlab("Vaping") +
#   ylab(paste0("VST-Expression (", DEResults_VapingLast6Mo$symbol[index_feature], ")")) +
#   scale_color_manual(values = palette_color_vaping) +
#   theme(legend.position = "none") +
#   facet_grid(~ RecruitmentCenter) + # facet_grid(~ Age <= 7)
#   NULL
# 




# 
# library(ggbeeswarm)
# DEResults_VapingLast6Mo_age_ethn$pvalue %>% hist
# quickplot_vaping <- function(index_feature) {
# tmp_exp <- assay(vstcounts_naive)[DEResults_VapingLast6Mo_age_ethn$row[index_feature], ]
# ggplot(data = metadata_naive,
#        aes(VapingLast6Mo, tmp_exp, color = VapingLast6Mo)) +
#   geom_quasirandom(width = 0.2, alpha = 0.7) +
#   theme_few() +
#   xlab("Vaping in Last 6 Months") +
#   ylab(paste0("VST-Expression (", DEResults_VapingLast6Mo_age_ethn$symbol[index_feature], ")")) +
#   scale_color_manual(name = "Vaping in Last\n6 Months", values = palette_color_vaping) +
#   # facet_grid( ~ RecruitmentCenter) +
#   theme(legend.position = "none")
# }
# 
# # 600 x 500
# plot_grid(quickplot_vaping(1), quickplot_vaping(2), quickplot_vaping(3), quickplot_vaping(4))
# 
# # 400 x 300
# quickplot_vaping(1) + facet_grid(~ RecruitmentCenter)
# quickplot_vaping(4) + facet_grid(~ RecruitmentCenter)
# 
# 
# quickplot_vaping(which(DEResults_VapingLast6Mo_age_ethn$symbol == "CCDC85A")) + facet_grid(~ RecruitmentCenter)
# quickplot_vaping(which(DEResults_VapingLast6Mo_age_ethn$symbol == "EVI2A")) + facet_grid(~ RecruitmentCenter)
# 
# 
# DEResults_VapingLast6Mo_center %>% head
# 
# 
# # ggplot(
# #   data = metadata_naive %>% filter(SID != "102"),
# #   aes(x = VapingPast30,
# #       y = assay(vstcounts)[DEResults_VapingPast30$row[6], ],
# #       color = VapingPast30)) +
# #   geom_quasirandom(alpha = 0.8, size = 2, width = 0.2) +
# #   theme_few()
# 
# ggplot(
#   data = metadata_naive %>% filter(SID != "102"),
#   aes(
#     x = Fres,
#     y = assay(vstcounts_naive)[DEResults_Fres$row[4], metadata_naive$SID != "102"],
#     color = VapingLast6Mo
#   )) +
#   geom_point(alpha = 0.8, size = 2) +
#   theme_few()
# 
# ggplot(
#   data = metadata_naive,
#   aes(
#     x = `FVC%`,
#     y = assay(vstcounts)[DEResults_FVCpp$row[1], ],
#     color = VapingLast6Mo
#   )
# ) +
#   geom_point(alpha = 0.8, size = 2) +
#   theme_few()
# 
# ggplot(
#   data = metadata_naive,
#   aes(
#     x = `Exhaled NO`,
#     y = assay(vstcounts)[DEResults_eNO$row[2], ],
#     color = VapingLast6Mo
#   )
# ) +
#   geom_point(alpha = 0.8, size = 2) +
#   theme_few()
# 
# 
# DEResults_VapingLast6Mo_age_ethn %>%
#   filter(padj < 0.1) %>%
#   View()
# 
# 
# 
# # 
# # object <- vstcounts
# # rv <- rowVars(assay(object)[1:100, ])
# # select <-
# #   order(rv, decreasing = TRUE)[seq_len(min(1, length(rv)))]
# # pca <- prcomp(t(assay(object)[select,]))
# 
# 
# filter_genes_variance <- 
#   vstcounts_naive %>%
#   assay %>%
#   apply(., 1, mean) %>%
#   order(decreasing = T) %>%
#   .[1:5000]
# pca_naive <-
#   assay(vstcounts_naive) %>%
#   .[filter_genes_variance, ] %>%
#   t %>%
#   prcomp(scale = T)
# percentexp <-
#   formatC(pca_naive$sdev^2 / sum(pca_naive$sdev^2) * 100, digits = 1, format = "f") %>%
#   paste0("PC", 1:length(.), " (", ., "%)")
# metadata_for_pca_naive <-
#   metadata_naive %>%# filter(SID != "102") %>%
#   bind_cols(PC1 = pca_naive$x[ , 1],
#             PC2 = pca_naive$x[ , 2])
# 
# 
# metadata_naive %>% filter(NewID %in% # odd QC metrics
#                             c("Sample18", "Sample12")) %>% .$SID
# legend_title <-
#   "Vaping in\nLast 6 Months"
# ggplot(metadata_for_pca_naive,
#        aes(PC1, PC2,
#            shape = VapingLast6Mo, color = VapingLast6Mo)) +
#   geom_point(size = 4, alpha = 0.75) +
#   geom_text_repel(aes(
#     label = SID),
#     # label = if_else(PC2 < -50, as.character(SID), "")),
#     show.legend = F) +
#   # facet_grid(~ RecruitmentCenter) +
#   theme_few() + theme(legend.position = "bottom") +
#   xlab(percentexp[1]) + ylab(percentexp[2]) +
#   scale_color_manual(legend_title, values = palette_color_vaping) +
#   scale_shape_manual(legend_title, values = palette_shape_vaping) +
#   NULL 
# legend_title <-
#   "Vaping in\nLast 6 Months"
# ggplot(metadata_for_pca,
#        aes(PC1, PC2,
#            shape = VapingLast6Mo, color = VapingLast6Mo)) +
#   geom_point(size = 4, alpha = 0.75) +
#   geom_text_repel(aes(label = if_else(PC2 < -50, as.character(SID), "")),
#                   show.legend = F) +
#   facet_grid(~ RecruitmentCenter) +
#   theme_few() + theme(legend.position = "bottom") +
#   xlab(percentexp[1]) + ylab(percentexp[2]) +
#   scale_color_manual(legend_title, values = palette_color_vaping) +
#   scale_shape_manual(legend_title, values = palette_shape_vaping) +
#   NULL 
# 
# 
# 
# 
# 
# 
# 
# left_join(DEResults_VapingLast6Mo_center, DEResults_VapingLast6Mo_cityaurora, by = "row") %>%
#   ggplot(., aes(log2FoldChange.x, log2FoldChange.y)) +
#   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point(alpha = 0.25)
# 
# 
# left_join(DEResults_VapingLast6Mo_center, DEResults_VapingLast6Mo_cityaurora, by = "row") %>%
#   ggplot(., aes(pvalue.x, pvalue.y)) +
#   geom_abline(slope = 1, intercept = 0) + 
#   geom_point(alpha = 0.5) +
#   scale_y_log10() + scale_x_log10()
# 
# 
# # 450 x 350
# ggplot(metadata_for_pca_SIDfiltered %>% filter(SID != "102"),
#        # mutate(`Race (Non-White)` = Race_white %>% as.logical() %>% `!`),
#        # mutate(`Latino (Ethnicity)` = Latino %>% `==`(., 2) %>% as.logical()),
#        aes(PC1, PC2,
#            # shape = VapingLast6Mo,
#            color = RecruitmentCenter)) +
#   geom_point(size = 4, alpha = 0.75) +
#   geom_text_repel(aes(label = if_else(PC2 < -50, as.character(SID), "")),
#                   show.legend = F) +
#   theme_few() + theme(legend.position = "bottom") +
#   xlab(percentexp[1]) + ylab(percentexp[2]) +
#   scale_color_brewer(palette = "Dark2") +
#   # scale_color_manual("Center", values = palette_color_center) +
#   # scale_shape_manual(legend_title, values = palette_shape_vaping) +
#   NULL
# 
# 
# # 450 x 350
# legend_title <-
#   "Vaping Last 30 Days (n = 5)"
# ggplot(metadata_for_pca_SIDfiltered %>% filter(SID != "102"),
#        aes(PC1, PC2,
#            shape = VapingPast30,
#            color = VapingPast30)) +
#   geom_point(size = 4, alpha = 0.75) +
#   geom_text_repel(aes(label = if_else(PC1 < -33, as.character(SID), "")),
#                   show.legend = F) +
#   theme_few() + theme(legend.position = "bottom") +
#   xlab(percentexp[1]) + ylab(percentexp[2]) +
#   scale_color_manual(legend_title, values = palette_color_vaping) +
#   scale_shape_manual(legend_title, values = palette_shape_vaping) +
#   NULL
# 
# 
# 
# 
# 
# 
# 
