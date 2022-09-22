This folder contains all the code.  

Details about the files in this folder:

File | Description
---|---------------------------------------------------------------------
01_sesame_getbeta_m | This is the file used to run the SeSAMe pipeline. Beta- and M-values from this file are written out to DataProcessed/methylation.
01a_sesame_getsex | This is the file used to complete the pre-processing check of biological sex inferred from methylation vs. clinical sex. The output is a data frame with inferred sex and all chromosomal intensities.
 01b_sesameQC | This file completes diagnostic checks for dye-bias and beta distributions for the SeSAMe pipeline.
01c_sesame_custom_pipeline | The purpose of this code is to run an alternative version of the "noob + BMIQ" procedure in the SeSAMe pipeline, but with the dye bias correction step moved to the end of the pipeline. Some diagnostic plots are created as well. 
 02_minfi_noobswan | This file was run in order to compare the pre-processing quality of the minfi pipeline (noob + 'SWAN' method) to the SeSAMe pipeline.
