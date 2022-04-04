library(edgeR)

final_counts <-
  counts_HC[filter_gene_cpm, ]

d <- DGEList(final_counts)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
disp <- d$tagwise.dispersion


mu <- apply(final_counts, 1, mean)

library(ssizeRNA)
powercalc <-
  ssizeRNA_vary(
  nGenes = nrow(final_counts), pi0 = 0.9, m = 50, mu = mu,
  disp = disp, fc = 1.5, fdr = 0.05,
  power = 0.8, maxN = 50)

50*0.9
powercalc$power
