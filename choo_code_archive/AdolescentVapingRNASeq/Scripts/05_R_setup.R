dyn.load("/home/liucu/bin/libiconv1.16/lib/libiconv.so")
dyn.load("/home/liucu/bin/libiconv1.16/lib/libiconv.so.2")

Sys.setenv(XML_CONFIG = "/home/liucu/libxml2-2.7.2/xml2-config")

export LD_LIBRARY_PATH=/home/liucu/libxml2-2.7.2/:/usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64


R
setwd("/home/biostats_share/liucu/ObeseAsthma_RNAseq/")

library(tidyverse)
library(DESeq2)