# libraries ====================================================================


library(tidyverse)
library(readxl)
library(choomisc)
library(skimr)
library(caret)
library(ggridges)
library(ggthemes)
library(writexl)
library(ggsignif)
library(cowplot)
library(ggbeeswarm)

library(tidyverse)
library(DESeq2)
library(data.table)

library(ggrepel)


library(RUVSeq)
library(vegan)




# palettes =====================================================================

palette_color_vaping <-
  c(Control = "black", Vaping = "#b2182b")
palette_color_center <-
  c(Pueblo = "#e41a1c", `CommCity/Denver` = "#377eb8", Aurora = "#4daf4a")
palette_shape_vaping <-
  c(Control = 1, Vaping = 19)
