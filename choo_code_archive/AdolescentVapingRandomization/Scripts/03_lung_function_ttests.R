palette_vape <-
  c(Vaping = "#C23151", Control = "#808080")
# vapinglabel <- "Vaping in\nLast 6 Months"
vapinglabel <- ""


lm(R5 ~ vapinggroup + Age + I(Gender == "1"), metadata) %>% summary
lm(R5 ~ vapinggroup + Age, metadata) %>% summary
t.test(R5 ~ vapinggroup, metadata)


# by vaping group
ggplot(metadata,
       aes(vapinggroup, R5, color = vapinggroup)) +
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
lm(R5 ~ vapinggroup + Age, metadata) %>% summary
plot_grid(
ggplot(metadata,
       aes(Age + 7, R5, color = vapinggroup)) +
  geom_point(size = 1.5, alpha = 0.8) +
  # geom_smooth(aes(group = 1), method = 'loess', se = F,
  #             size = 0.5, color = "#000000") +
  stat_summary(geom = "errorbar", width = 0.25, color = "black") +
  stat_summary(fun = mean, fun.max = mean, geom = "errorbar",
               width = 0.5, color = "black") +
  theme_few() +
  xlab("Age (years)") +
  scale_color_manual(name = vapinglabel, values = palette_vape) +
  theme(legend.position = "none"),
# by gender
ggplot(metadata,
       aes(factor(Gender), R5, color = vapinggroup)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(name = vapinglabel, values = palette_vape) +
  scale_x_discrete("Gender", labels = c("M", "F", "Other")) +
  stat_summary(geom = "errorbar", width = 0.25, color = "black") +
  stat_summary(fun = mean, fun.max = mean, geom = "errorbar",
               width = 0.5, color = "black") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_few(),
rel_widths = c(4, 5)
)

ggsave(filename = "20210113_cov.png", width = 8, height = 3)
system('open 20210113_cov.png')






variables_to_test <- c("R5", "X5", "VT", "AX", "Fres")
plots_to_export <-
  variables_to_test %>%
  sort %>%
  lapply(., function(feature) {

        pval <- t.test(as.formula(paste0(feature, " ~ vapinggroup")), metadata)$p.value %>%
          formatC(digits = 2, flag = "#") %>% 
          paste0("P = ", .)
        
      ggplot(metadata,
           aes_string(x = "vapinggroup", y = feature, color = "vapinggroup")) +
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
  plots_to_export[[4]]
  )
# ggsave(filename = "20210113_lungfxn.png", width = 5, height = 5)
# system('open 20210113_lungfxn.png')



sapply(variables_to_test,
       function(feature) {
         lm(as.formula(paste0(feature, " ~ vapinggroup + Age + I(Gender == 1)")),
            metadata) %>% summary %>% coef %>% .[2, 4] %>%
           formatC(digits = 2, flag = "#") %>% 
           paste0("P = ", .)
       })
