palette_vape <-
  c(Vaping = "#C23151", Control = "#808080")
# vapinglabel <- "Vaping in\nLast 6 Months"
vapinglabel <- ""


lm((X20 )~ VapingLast6Mo + Age + I(Gender == "1"), metadata) %>% summary
lm(R5 ~ VapingLast6Mo + Age, metadata) %>% summary
t.test(R5 ~ VapingLast6Mo, metadata)


# by vaping group
ggplot(metadata_all,
       aes(VapingLast6Mo, R5, color = VapingLast6Mo)) +
  geom_point() +
  scale_color_manual(name = vapinglabel, values = palette_vape) +
  xlab(vapinglabel) +
  stat_summary(geom = "errorbar", width = 0.25, color = "black") +
  stat_summary(fun = mean, fun.max = mean, geom = "errorbar",
               width = 0.5, color = "black") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_few() + theme(legend.position = "none") +
  geom_signif(annotations = "P = 0.033", comparisons = list(c("Control", "Vaping")))

# by age
lm(R5 ~ VapingLast6Mo + Age, metadata) %>% summary

# "#1c9099"

plot_grid(
  ggplot(metadata_all,
         aes(Age + 7, VT)) +
    geom_point(size = 1.5, alpha = 0.8, color = "#d95f0e") +
    # geom_smooth(aes(group = 1), method = 'loess', se = F,
    #             size = 0.5, color = "#000000") +
    stat_summary(geom = "errorbar", width = 0.25, color = "black") +
    stat_summary(fun = mean, fun.max = mean, geom = "errorbar",
                 width = 0.5, color = "black") +
    theme_few() +
    xlab("Age (years)") +
    # scale_color_manual(name = vapinglabel, values = palette_vape) +
    theme(legend.position = "none"),
  # by gender
  ggplot(metadata_all,
         aes(factor(RecruitmentCenter), VT)) +
    geom_point(size = 1.5, alpha = 0.8, color = "#d95f0e") +
    # scale_color_manual(name = vapinglabel, values = palette_vape) +
    scale_x_discrete("Recruitment Center", labels = c("A", "C/D", "P")) +
    stat_summary(geom = "errorbar", width = 0.25, color = "black") +
    stat_summary(fun = mean, fun.max = mean, geom = "errorbar",
                 width = 0.5, color = "black") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme_few() + theme(legend.position = "none"),
  # by gender
  ggplot(metadata_all,
         aes(factor(Gender), VT)) +
    geom_point(size = 1.5, alpha = 0.8, color = "#d95f0e") +
    # scale_color_manual(name = vapinglabel, values = palette_vape) +
    scale_x_discrete("Gender", labels = c("M", "F", "NB")) +
    stat_summary(geom = "errorbar", width = 0.25, color = "black") +
    stat_summary(fun = mean, fun.max = mean, geom = "errorbar",
                 width = 0.5, color = "black") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme_few(),
  rel_widths = c(4, 4, 4),
  nrow = 1
)
# 600 x 250

# ggsave(filename = "20210113_cov.png", width = 8, height = 3)
# system('open 20210113_cov.png')







variables_to_test <- c("R5", "R20", "R20 - R5", "X5", "X20", "X20 - X5", "VT", "AX", "Fres")
plots_to_export <-
  variables_to_test %>%
  sort %>%
  lapply(., function(feature) {
    
    pval <-
      # wilcox.test(as.formula(paste0(feature, " ~ VapingLast6Mo")), metadata_all)$p.value %>%
      # t.test(as.formula(paste0(feature, " ~ VapingLast6Mo")), metadata_all)$p.value %>%
      lm(as.formula(paste0(feature, " ~ VapingLast6Mo + RecruitmentCenter")),
         metadata_all,# %>% filter(RecruitmentCenter == "Pueblo")
         ) %>% summary %>% coef %>% .[2, 4] %>%
      formatC(digits = 2, flag = "#") %>% 
      paste0("P = ", .)
    
    ggplot(metadata_all,# %>% filter(RecruitmentCenter == "Pueblo"),
           aes_string(x = "VapingLast6Mo", y = feature, color = "VapingLast6Mo")) +
      geom_quasirandom(width = 0.2) +
      scale_color_manual(name = vapinglabel, values = palette_vape) +
      xlab(vapinglabel) +
      stat_summary(#fun.ymin = function(x) { quantile(x, 0.25) },
        #fun.ymax = function(x) { quantile(x, 0.75) },
        geom = "errorbar", width = 0.25, color = "black") +
      stat_summary(fun = mean, fun.max = mean, geom = "errorbar",
                   width = 0.5, color = "black") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                         expand = c(0.15, 0)) +
      theme_few() + theme(legend.position = "none") +
      geom_signif(annotations = pval, margin_top = 0.15, textsize = 3,
                  comparisons = list(c("Control", "Vaping")))
    
    
  })

plot_grid(
  plots_to_export[[1]],
  plots_to_export[[2]],
  plots_to_export[[3]],
  plots_to_export[[4]],
  plots_to_export[[5]],
  plots_to_export[[6]],
  plots_to_export[[7]],
  plots_to_export[[8]],
  plots_to_export[[9]]
)



plot_grid(
  plots_to_export_NOSIGNIF[[1]],
  plots_to_export_NOSIGNIF[[2]],
  plots_to_export_NOSIGNIF[[3]],
  plots_to_export_NOSIGNIF[[4]],
  plots_to_export_NOSIGNIF[[5]],
  plots_to_export_NOSIGNIF[[6]],
  plots_to_export_NOSIGNIF[[7]],
  plots_to_export_NOSIGNIF[[8]],
  plots_to_export_NOSIGNIF[[9]]
)

extract_pvalue_IOS <-
  function(feature, rhs) {
    lm(as.formula(paste0(feature, " ~ ", rhs)),
     metadata_all) %>% summary %>% coef %>% .[2, 4]
  }

IOS_pvalue_table <-
  variables_to_test %>%
  sapply(function(feature) {
    c(Age = extract_pvalue_IOS(feature, "Age"),
      Gender = extract_pvalue_IOS(feature, "GenderMale"),
      Ethnicity = extract_pvalue_IOS(feature, "EthnicityLatino"),
      Center = 
        anova(lm(as.formula(paste0(feature, " ~ RecruitmentCenter")), metadata_all),
              lm(as.formula(paste0(feature, " ~ 1")), metadata_all)) %>% .$`Pr(>F)` %>% .[2],
    Vaping_NoAdj = extract_pvalue_IOS(feature, "VapingLast6Mo"),
    Vaping_AdjCenter = extract_pvalue_IOS(feature, "VapingLast6Mo + RecruitmentCenter"),
    Vaping_AdjAgeEthnGender = extract_pvalue_IOS(feature, "VapingLast6Mo + Age + EthnicityLatino + GenderMale"),
    Vaping_PuebloOnly = lm(as.formula(paste0(feature, " ~ VapingLast6Mo")), metadata_all %>% filter(RecruitmentCenter == "Pueblo")) %>%
      summary %>% coef %>% .[2, 4])
      })

IOS_pvalue_table %>% t %>%
  as_tibble(rownames = "Feature") %>% mutate_at(.vars = 2:ncol(.), .funs = function(x) { if_else(x < 0.1, paste0(formatC(x, digits = 2), "*"), formatC(x, digits = 2))}) %>%
  View

plots_to_export_NOSIGNIF <-
  variables_to_test %>%
  sort %>%
  lapply(., function(feature) {
    
    pval <-
      # wilcox.test(as.formula(paste0(feature, " ~ VapingLast6Mo")), metadata)$p.value %>%
      # t.test(as.formula(paste0(feature, " ~ VapingLast6Mo")), metadata)$p.value %>%
      lm(as.formula(paste0(feature, " ~ VapingLast6Mo + RecruitmentCenter")),
         metadata# %>% filter(RecruitmentCenter == "Pueblo")
      ) %>% summary %>% coef %>% .[2, 4] %>%
      formatC(digits = 2, flag = "#") %>% 
      paste0("P = ", .)
    
    ggplot(metadata_all,# %>% filter(RecruitmentCenter == "Pueblo"),
           aes_string(x = "VapingLast6Mo", y = feature, color = "VapingLast6Mo")) +
      geom_quasirandom(width = 0.2) +
      scale_color_manual(name = vapinglabel, values = palette_vape) +
      xlab(vapinglabel) +
      stat_summary(#fun.ymin = function(x) { quantile(x, 0.25) },
        #fun.ymax = function(x) { quantile(x, 0.75) },
        geom = "errorbar", width = 0.25, color = "black") +
      stat_summary(fun = mean, fun.max = mean, geom = "errorbar",
                   width = 0.5, color = "black") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                         expand = c(0.15, 0)) +
      theme_few() + theme(legend.position = "none") +
      geom_signif(annotations = pval, margin_top = 0.15, textsize = 3,
                  comparisons = list(c("Control", "Vaping")))
    
    
  })
plots_to_export_NOSIGNIF[[5]] + facet_grid(~ RecruitmentCenter) 
plots_to_export[[7]] + facet_grid(~ RecruitmentCenter) 
plots_to_export_NOSIGNIF[[5]] + facet_grid(~ RecruitmentCenter) 
ggsave(filename = "20210326_lungfxn.png", width = 8, height = 8)
system('open 20210326_lungfxn.png')



sapply(variables_to_test,
       function(feature) {
         lm(as.formula(paste0(feature, " ~ VapingLast6Mo + Age + I(Gender == 1)")),
            metadata_all) %>% summary %>% coef %>% .[2, 4] %>%
           formatC(digits = 2, flag = "#") %>% 
           paste0("P = ", .)
       })
