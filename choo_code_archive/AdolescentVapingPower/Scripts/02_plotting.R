library(tidyverse)
library(ggthemes)

pwrresults <-
  tibble(percent =
           rep(percent_reduction_values * 100, 2),
         cohen = c(2.27*percent_reduction_values/sd_values_small,
                   2.27*percent_reduction_values/sd_values_big),
         N = c(n_per_group_small, n_per_group_big),
         SD_assume = rep(c("small", "large"), each = length(percent_reduction_values))
         ) %>%
  group_by(SD_assume) %>%
  group_split() %>%
  map(~ mutate(.x,
               N_label = if_else(!duplicated(ceiling(N)), as.character(ceiling(N)), ""))) %>%
  bind_rows()

ggplot(
  data = pwrresults,
  aes(x = percent, y = N, color = SD_assume)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = N_label), show.legend = F) +
  geom_line() +
  theme_few() +
  xlab("% Reduction") +
  ylab("N/Arm") +
  scale_color_brewer(name = "Assumed SD", palette = "Set2") +
  theme(legend.position = c(0.8, 0.8))

ggsave(
  'powerfig.png', width = 4, height = 3, dpi = 200
)

system('open powerfig.png')
