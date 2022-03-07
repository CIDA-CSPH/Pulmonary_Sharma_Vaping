library(ComplexHeatmap)

filter_samples <-
  (metadata_final$SID != "102") 

metadata_final_filter <-
  metadata_final[filter_samples, ]

vstcounts <-
  DESeqDataSetFromMatrix(
    countData = final_counts_mat[ , filter_samples],
    colData = metadata_final_filter,
    design = ~ VapingLast6Mo) %>%
  vst()
vst_dm <-
  vstcounts %>%
  assay(.) %>%
  .[filter_genes_variance, ] %>% t %>% dist
library(vegan)
adonis(vst_dm ~ Age + I(Gender == "1") + I(Latino == "1") + RecruitmentCenter, #+ VapingLast6Mo,
       metadata_final_filter) %>%
  .$aov.tab %>%
  View()
adonis(vst_dm ~ Age + Latino + VapingLast6Mo,
       metadata_final_filter) %>%
  .$aov.tab %>%
  View()



resids_no_group <-
  deseqobj_VapingLast6Mo_RUV %>%
  vst() %>%
  assay %>%
  t %>%
  # quantile_normalisation %>%
  apply(., 2,
        function(y) {
          lm(y ~ AgeScaled + EthnicityLatino + GenderMale +
               RUV1 + RUV2 + RUV3 + RUV4,
             metadata_final) %>% resid
        })



# heatmap_mat <-
#   deseqobj_VapingLast6Mo_RUV %>%
#   vst() %>%
#   assay  %>%
#   .[
#     DEResults_VapingLast6Mo_RUV %>% 
#       # filter(padj < 0.01) %>%
#       filter(padj < 0.05 & abs(log2FoldChange > sqrt(2))) %>%
#       .$row,
#   ] %>%
#   apply(., 1, function(x) { (x - mean(x))/sd(x)}) %>%
#   # apply(., 1, function(x) { (x - min(x))/(max(x) - min(x))}) %>%
#   `rownames<-`(NULL) %>% `colnames<-`(NULL)


heatmap_mat <- 
  resids_no_group %>% 
  .[ ,
    DEResults_VapingLast6Mo_RUV %>% 
      # filter(padj < 0.01) %>%
      filter(padj < 0.05 & abs(log2FoldChange) > sqrt(2)) %>%
      .$row
  ] %>%
  # t() %>%
  quantile_normalisation %>%
  apply(., 2, function(x) { (x - min(x))/(max(x) - min(x))}) %>%
  # apply(., 1, function(x) { (x - mean(x))/sd(x)}) %>%
  `rownames<-`(NULL) %>%
  `colnames<-`(NULL)

# https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}


library(viridis)
Heatmap(heatmap_mat,
        name = "Relative\nExpression",
        col = viridis(200, option = "B"),
        # clustering_distance_rows = "pearson",
        # clustering_distance_columns = "pearson",
        right_annotation =
          rowAnnotation(df = metadata_final %>% # dplyr::slice(-35) %>%
                          select(VapingLast6Mo) %>% as.data.frame(),
                        col = list(VapingLast6Mo = palette_color_vaping),
                        show_annotation_name = F)
        )





DEResults_VapingLast6Mo_RUV %>%
  mutate(dir = case_when(
    padj < 0.05 & (log2FoldChange > 0) ~ "up",
    padj < 0.05 & (log2FoldChange < 0) ~ "down",
    T ~ "N.S."
  )) %>%
  # filter(padj < 0.05) %>%
  ggplot(data = .,
         aes(log2FoldChange, -log10(pvalue),
             color = dir, shape = padj < 0.05)) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_point(alpha = 0.2) +
  geom_text_repel(aes(label = if_else(
    abs(log2FoldChange) > 3 | (pvalue < 1e-12 & log2FoldChange > 0) |
      (pvalue < 1e-25 & log2FoldChange < 0.5), symbol, "")),
    size = 3, segment.size = 0.2, box.padding = 0.1) +
  scale_color_manual(values = c("#1a9850", "black", "#3288bd")) +
  scale_shape_manual(values = c(1, 19)) +
  theme_few() + theme(legend.position = "none") +
  xlab(expression(log[2]~"(Fold Change)")) +
  ylab(expression(log[10]~"(P-Value)"))
