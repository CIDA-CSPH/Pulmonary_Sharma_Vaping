library(WebGestaltR)
dir.create("./Output/20210412_webgestalt")

ls(pattern = "DEResults_*")

export_rnk_for_gsea <-
  function(deresults_name) {
    get(deresults_name) %>%
      filter(!is.na(pvalue)) %>% select(symbol, log2FoldChange) %>%
      write_tsv(
        paste0("./Output/20210412_webgestalt/", deresults_name, ".rnk"),
        col_names = F)
}


export_rnk_for_gsea("DEResults_VapingLast6Mo_RUV")


# 
# webgest_wrapper <- 
#   function(filein, refin, outname) {
#     WebGestaltR(
#       interestGeneFile = filein,
#       interestGeneType = "genesymbol",
#       referenceGeneType = "genesymbol",
#       referenceGeneFile = refin,
#       enrichDatabase = c("geneontology_Biological_Process", "pathway_Panther", "pathway_Reactome"),
#       outputDirectory = "./Output/20210311_webgestalt", projectName = outname,
#       fdrThr = 1, reportNum = 20
#     )
#   }
# 
# webgest_wrapper("./Output/20210311_webgestalt/VapingLast6Mo_AE_target.txt",
#                 "./Output/20210311_webgestalt/VapingLast6Mo_AE_ref.txt",
#                 "VapingLast6Mo_AE")
# 
# 
# # GSEA
# 
# write_xlsx(
#   list(VapingLast6Mo = DEResults_VapingLast6Mo_age_ethn,
#        VapingLast30 = DEResults_VapingPast30_age_ethn,
#        VapingTimesDay = DEResults_VapingTimesPer_age_ethn,
#        LungFxnAX = DEResults_AX_age_ethn,
#        LungFxnFres = DEResults_Fres_age_ethn,
#        LungFxnR5 = DEResults_R5_age_ethn,
#        LungFxnX5 = DEResults_X5_age_ethn,
#        LungFxnVT = DEResults_VT_age_ethn) %>%
#     map(~ filter(.x, padj < 0.1)),
#   path = "./Output/20200311_Vaping_DE_ageethnicity.xlsx"
# )
# 
# write_xlsx(
#   list(
#   GO_AgeEthn = read_tsv("./Output/20210311_webgestalt/webresults/Project_wg_Vaping_AE_GOBP/enrichment_results_wg_result1615509892.txt"),
#   Reactome_AgeEthn = read_tsv("./Output/20210311_webgestalt/webresults/Project_wg_Vaping_AE_Reactome/enrichment_results_wg_result1615507746.txt"),
#   GO_Center = read_tsv("./Output/20210311_webgestalt/webresults/Project_wg_Vaping_Center_GOBP/enrichment_results_wg_result1615509560.txt"),
#   Reactome_Center = read_tsv("./Output/20210311_webgestalt/webresults/Project_wg_Vaping_Center_Reactome/enrichment_results_wg_result1615508674.txt")),
#   path = "./Output/20200311_Vaping_GSEAResults.xlsx")
# 
# system("open ./Output/20200311_Vaping_GSEAResults.xlsx")
# 
# 






GSEAResults_Vaping <-
  read_tsv("./Output/20210412_webgestalt/Project_wg_result1618254369_Reactome/enrichment_results_wg_result1618254369.txt") %>%
  transmute(ReactomeID = geneSet,
            Description = description,
            NES = normalizedEnrichmentScore,
            `p-value` = pValue,
            FDR,
            `Leading Edge` = userId)
