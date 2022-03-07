library(pwr)

# est population mean = 2.27, sd = 3.84
run_pwr_calc <-
  function (percent_reduction, sd = 3.84) {

pwr.t.test(
  d = (2.27 * percent_reduction)/sd,
  power = 0.8,
  sig.level = 0.05,
  type = "paired", # "one.sample"
  alternative = "two.sided")

}

percent_reduction_values <-
  seq(0.15, 1, 0.05)
set.seed(1234)

sd_values_small <-
  sd(2.27 - runif(n = 150, -2, 2)/4)
sd_values_big <-
  sd(2.27 - runif(n = 150, -4, 4)/4)

2.27*percent_reduction_values
percent_reduction_values/sd_values

n_per_group_small <-
  sapply(percent_reduction_values,
       function(x) { run_pwr_calc(percent_reduction = x,
                                  sd = sd_values_small)$n })

n_per_group_big <-
  sapply(percent_reduction_values,
         function(x) { run_pwr_calc(percent_reduction = x,
                                    sd = sd_values_big)$n })
