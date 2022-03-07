library(limma)

exprs_obj <-
  methyl_mvals[filter_probes, -1] %>%
  as.matrix
probes_tested <-
  methyl_mvals$Probe[filter_probes]
limma_obj <- 
  lmFit(exprs_obj, 
        design = model.matrix(~ 1 + VapingLast6Mo + Age + GenderMale + RecruitmentCenter + EthnicityLatino, final_metadata))
limma_results <-
  eBayes(limma_obj) %>% limma::topTable(limma_results, coef = 2, number = 1e6) %>%
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
bacon_adj <- bacon::bacon(teststatistics = limma_results$t.1)
bacon_adj@teststatistics %>% head
bacon_P <- pval(bacon_adj)[ , 1]
bacon_FDR <- p.adjust(bacon_P, "fdr")

# plot(bacon_adj, type = "qq")

limma_results %>%
  as_tibble(rownames = "probe_index") %>% 
  transmute(probe_index = probe_index,
            probe = as.numeric(probe_index) %>% probes_tested[.],
         Effect = coefficients.VapingLast6MoVaping) %>%
  bind_cols(P = bacon_P, FDR = bacon_FDR) %>%
  left_join(., EPIC.hg38.manifest_df %>% select(probe, gene), by = "probe") %>%
  arrange(P)

ggplot(data = final_metadata,
       aes(x = VapingLast6Mo, y = exprs_obj[2804, ])) +
  geom_point()

