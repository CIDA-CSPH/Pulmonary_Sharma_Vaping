R
library(tidyverse)
library(data.table)

setwd('/home/biostats_share/liucu/Vaping/')

star_filepaths <-
  list.files("aligned_reads", pattern = "*ReadsPerGene.out.tab",
             recursive = T, full.names = T)
sample_names <-
  star_filepaths %>%
  str_split(., pattern = "_") %>%
  map_chr(~ .[3]) %>%
  str_pad(string = ., width = 2, pad = "0", side = "left") %>%
  paste0("Sample", .)

gene_names <-
  read_tsv(star_filepaths[1], col_names = F) %>% select(1) %>% .$X1

counts <-
  map(star_filepaths, function(x) { fread(x) %>% .[ , V3] }) %>% # V3 = sense strand
  set_names(sample_names) %>%
  as_tibble %>%
  select(names(.) %>% sort) %>%
  bind_cols(Feature = gene_names, .) %>%
  slice(-(1:4)) %>%
  arrange(Feature)


dir.create("Output")
write_tsv(counts, "./Output/20210307_counts.txt")