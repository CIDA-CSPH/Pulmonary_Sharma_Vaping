This folder contains all the code.  

Details about the files in this folder:

File | Description
---|---------------------------------------------------------------------
01_sesame_getbetas | This is the file used to run the SeSAMe pipeline. Beta- and M-values from this file are written out to DataProcessed/methylation.
01a_sesame_getsex | This is the file used to complete the pre-processing check of biological sex inferred from methylation vs. clinical sex. The output is a data frame with inferred sex and all chromosomal intensities.
 01b_sesameQC | This file completes diagnostic checks for dye-bias and beta distributions for the SeSAMe pipeline as well as MDS plots to look at sample clustering.
01c_sesame_custom_pipeline | The purpose of this code is to run an alternative version of the "noob + BMIQ" procedure in the SeSAMe pipeline, but with the dye bias correction step moved to the end of the pipeline. Some diagnostic plots are created as well.
 02_RUVm | Run this file in order to run Removal of Unwanted Variability for Methylation (RUVm) as well as look at sample clustering after RUVm.
 02a_RUVr | This file attempts to use RUVr, but for the methylation data in this report. 
