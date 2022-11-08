# Check Correlations ------------------------------------------------------

check_corr<- function(y, x, dat){
  #T.Test for categorical with 2 levels
  if(is.factor(dat[[x]]) & length(levels(dat[[x]])) == 2){
    temp_res <- t.test(dat[[y]] ~ dat[[x]])
    pval <- temp_res$p.value
  }
  
  #One-way AOV for categorical with > 2 levels
  if(is.factor(dat[[x]]) & length(levels(dat[[x]])) > 2){
    temp_res <- summary(aov(dat[[y]] ~ dat[[x]]))
    pval <- as.numeric(temp_res[[1]]$`Pr(>F)`[1])
  }
  
  #Two continuous variables
  if(is.numeric(dat[[x]])){
    temp_res <- cor.test(dat[[x]], dat[[y]])
    pval <- temp_res$p.value
  }
  return(pval)  
}

#Variables of interest
vars_of_interest <- c("recruitment_center", "vape_status", "age")
ruv_k1_corr <- sapply(vars_of_interest, function(x) check_corr("ruv_k1", x, dat = clin_metadata))
ruv_k2_corr <- sapply(vars_of_interest, function(x) check_corr("ruv_k2", x, dat = clin_metadata))

cbind(ruv_k1_corr, ruv_k2_corr) %>% 
  as.data.frame() %>% 
  mutate(test_type = c("One-way AOV", "T-test", "Pearson Correlation")) %>% 
  dplyr::rename('ruv_K1' = ruv_k1_corr,
                "ruv_K2" = ruv_k2_corr,
                "Test Type" = test_type) %>% 
  gt::gt(rownames_to_stub = T) %>% 
  gt::fmt_number(columns = c(ruv_K1, ruv_K2), decimals = 3) %>% 
  gt::tab_style(style = gt::cell_text(color = "red"),
                locations = list(gt::cells_body(columns = "ruv_K2", rows = ruv_K2 < 0.05))
  ) %>% 
  gt::tab_options(table.width = 700) %>% 
  gt::cols_align(align = "left") %>% 
  gt::tab_header(
    title = "P-Values for Association"
  )
