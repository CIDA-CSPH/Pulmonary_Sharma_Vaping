library(RUVSeq)

filter_probes_anyNA <-
  methyl_mvals %>% apply(., 1, function(x) { is.na(x) %>% any})

filter_probes_autosomal <-
  EPIC.hg38.manifest_df$seqnames %in% paste0("chr", 1:22)

# filter_beta_range <- 
#   probes_beta_range > 0.1

filter_beta_range <- 
  order(probes_beta_range) %in% c(1:5000)

filter_probes <- !filter_probes_anyNA & filter_probes_autosomal & filter_beta_range


exprs_obj <-
  methyl_mvals[filter_probes, -1] %>%
  as.matrix
probes_tested <-
  methyl_mvals$Probe[filter_probes]

preruv_residuals <- 
  apply(exprs_obj, 1,
        function(y) { lm(y ~ EthnicityLatino + Age + GenderMale + VapingLast6Mo,
                         final_metadata) %>% 
            residuals()})


dm_preruvresid <- 
  preruv_residuals %>% dist(method = "euclidean")


plot(1:10,
     sapply(1:10,
            function(kval) {
ruv_results <-
  RUVr(
  x = exprs_obj,
  k = kval,
  residuals = t(preruv_residuals),
  isLog = T)

r2 <-
  adonis(dm_preruvresid ~ ruv_results$W)
  r2$aov.tab$R2[1]
}))



ruv_results <-
  RUVr(
    x = exprs_obj,
    k = 4,
    residuals = t(preruv_residuals),
    isLog = T)






pca_results <-
  prcomp(ruv_results$normalizedCounts %>% t, center = T, scale = F)
pca_percent_exp <-
  (pca_results$sdev^2) %>% `/`(., sum(.)) %>% `*`(100) %>% formatC(digits = 2) %>%
  gsub(" ", "", .) %>%
  paste0(" (", ., "%)") %>%
  paste0("PC", 1:length(.), .)



pca_df <-
  final_metadata %>%
  bind_cols(pca_results$x[ , 1:3] %>% as_tibble)

ggplot(pca_df, aes(PC1, PC2, color = VapingLast6Mo)) + #, color = sex, shape = sex)) +
  geom_point(alpha = 0.75) +
  geom_text_repel(
    aes(label = SID),
    color = "black", size = 2) +
  # scale_shape_manual(name = "Sex", values = palette_shape_sex, labels = c("Female", "Male")) +
  # scale_color_manual(name = "Sex", values = palette_color_sex, labels = c("Female", "Male")) +
  xlab(pca_percent_exp[1]) + ylab(pca_percent_exp[2]) +
  ggtitle("Methylation PCA") +
  theme_few() +
  # facet_grid(~ RecruitmentCenter) +
  geom_blank()






