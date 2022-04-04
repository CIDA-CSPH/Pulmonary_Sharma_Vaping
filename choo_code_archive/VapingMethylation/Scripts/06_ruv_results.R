library(limma)

filter_probes_range_full <-
  probes_beta_range > 0.05
filter_probes <- !filter_probes_anyNA & filter_probes_autosomal & filter_probes_range_full


exprs_obj <-
  methyl_mvals[filter_probes, -1] %>%
  as.matrix
probes_tested <-
  methyl_mvals$Probe[filter_probes]
limma_obj <- 
  lmFit(exprs_obj, 
        design = model.matrix(~ 1 + VapingLast6Mo + Age + GenderMale + EthnicityLatino + ruv_results$W,
                              final_metadata))
limma_results <-
  eBayes(limma_obj) %>%
  limma::topTable(., coef = 2, number = 1e6) %>%
  as.data.frame()
top10 <-
  limma_results %>%
  row.names() %>%
  as.numeric %>%
  probes_tested[.]

EPIC.hg38.manifest_df %>%
  # filter(probe == "cg21570723") %>%
  filter(probe %in% top10) %>%
  .$gene_HGNC



library(bacon)
bacon_adj <- bacon::bacon(teststatistics = limma_results$t)
bacon_adj@teststatistics %>% head
bacon_P <- pval(bacon_adj)[ , 1]
bacon_FDR <- p.adjust(bacon_P, "fdr")

# plot(bacon_adj, type = "qq")

limma_results_table <-
  limma_results %>%
  as_tibble(rownames = "probe_index") %>% 
  transmute(probe_index = probe_index,
            probe = as.numeric(probe_index) %>% probes_tested[.],
            # Effect = coefficients.VapingLast6MoVaping) %>%
            Effect = logFC) %>%
  bind_cols(P = bacon_P, FDR = bacon_FDR) %>%
  left_join(., EPIC.hg38.manifest_df %>% select(probe, gene), by = "probe") %>%
  arrange(P)

postruv_resid <-
  apply(exprs_obj, 1,
        function(y) { lm(y ~ EthnicityLatino + Age + GenderMale + ruv_results$W,
                         final_metadata) %>% 
            residuals()})

ggplot(data = final_metadata,
       aes(x = VapingLast6Mo,
           # y = exprs_obj[ limma_results_table$probe_index[1] %>% as.numeric(), ])) +
           y = postruv_resid[ limma_results_table$probe_index[1] %>% as.numeric(), ])) +
  geom_point()


sum(limma_results_table$FDR < 0.05)

