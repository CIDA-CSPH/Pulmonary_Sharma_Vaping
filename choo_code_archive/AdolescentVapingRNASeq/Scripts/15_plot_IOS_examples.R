
ls(pattern = "DEResults*")



signif_in_both <-
  which((DEResults_VT_RUV %>% filter(padj < 0.05) %>% .$row)
          %in% (DEResults_VapingLast6Mo_RUV %>% filter(padj < 0.05) %>% .$row))

plot_gene_by_IOS_and_vaping <-
  function(index) {
    ggplot(data = metadata_final %>% filter(!is.na(VT)) %>% .[deseqobj_VT_RUV %>% assays %>% .[["weights"]] %>% .[signif_in_both[index], ] %>% `==`(., 1), ],
           aes(VT,
               # vstcounts_preruv[signif_in_both[index],
               #                  deseqobj_VT_RUV %>% assays %>% .[["weights"]] %>% .[signif_in_both[index], ] %>% `==`(., 1) %>% .[.] %>% names],
               ruvresults_VT_RUV$normalizedCounts[signif_in_both[index],
                                                  deseqobj_VT_RUV %>% assays %>% .[["weights"]] %>% .[signif_in_both[index], ] %>% `==`(., 1) %>% .[.] %>% names],
               color = VapingLast6Mo,
               shape = VapingLast6Mo)) +
      geom_point(size = 2)  +
      theme_few() +
      scale_color_manual(name = "Exposure", values = palette_color_vaping) +
      scale_shape_manual(name = "Exposure", values = palette_shape_vaping) +
      ylab(paste0("RUV-Expression\n(", DEResults_VT_RUV$symbol[signif_in_both[index]], ", log2FC = ",
                  DEResults_VT_RUV$log2FoldChange[signif_in_both[index]] %>% formatC(digits = 2), ")"))
    
  }

# DEResults_VT_RUV %>% head
# plot_grid(plot_gene_by_IOS_and_vaping(1) + geom_smooth(aes(group = 1), se = F),
#           plot_gene_by_IOS_and_vaping(2) + geom_smooth(aes(group = 1), se = F),
#           plot_gene_by_IOS_and_vaping(3) + geom_smooth(aes(group = 1), se = F),
#           plot_gene_by_IOS_and_vaping(4) + geom_smooth(aes(group = 1), se = F))

# 800 x 300
plot_grid(plot_gene_by_IOS_and_vaping(1),
          plot_gene_by_IOS_and_vaping(2))
# plot_gene_by_IOS_and_vaping(3) + facet_grid( ~ EthnicityLatino)
