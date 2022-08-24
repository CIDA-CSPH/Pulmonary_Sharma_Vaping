####################Load Libraries########################
library(tidyverse)
library(here)
library(kableExtra)
library(sesame) #ver 1.15.6
library(readxl)

#################### Read in IDAT ########################
folder_raw_dat <- here("DataRaw/methylation/RawIdat/")

methyl_raw <- lapply(searchIDATprefixes(folder_raw_dat),
                          readIDATpair)
################## Run openSesame pipeline step-by-step #################
methyl_ses_process <- lapply(methyl_raw, 
                               function(pfx) {prepSesame(pfx, "QCDPB")})
methyl_process_RGSet <- SigSetsToRGChannelSet(methyl_ses_process)
