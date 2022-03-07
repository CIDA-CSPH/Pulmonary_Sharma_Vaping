################################################### -
## Title: Power Analysis
## Author: Ryan Peterson
## PI: Sharma
## Date Created: Tue Jun 29 10:18:37 2021
################################################### -

library(tidyverse)
library(here)
library(CIDAtools)

### Aim 1 #### 
n1 <- 200*.85; n2 <- 100*.85

## SDs:

# 1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1884751/
# 2: https://www.journalpulmonology.org/en-impulse-oscillometry-spirometry-passive-smoking-articulo-S2173511517300982
# 3: https://ejb.springeropen.com/articles/10.1186/s43168-020-00037-8
sd <- c(.13, .2, .25)
power = c(.7,.8,.9,.95)
d <- c(.2, .5, .8)
sig.level <- c(.05, .01, .005)

tab1 <- vec_power(pwr::pwr.t2n.test, n1 = n1, n2=n2, power = power) 

tab1 %>% 
  mutate(sd1 = d*sd[1], sd2 = d*sd[2], sd3 = d*sd[3])

vec_power(pwr::pwr.t2n.test, n1 = n1, n2=n2, d= d) %>% 
  mutate(sd1 = d*sd[1], sd2 = d*sd[2], sd3 = d*sd[3])





### Aim 3 ####

# Effect size	d	Reference
# Very small	0.01	[9]
# Small	0.20	[8]
# Medium	0.50	[8]
# Large	0.80	[8]
# Very large	1.20	[9]
# Huge	2.0	[9]

dropout <- .15

pwr::pwr.t.test(n = floor(75 * (1-dropout)), power= .8)

pwr::pwr.t.test(n = floor(75 * (1-dropout)), d = .2)

power.t.test(n = floor(75 * (1-dropout)), delta = 11.45, sd = 22.9)
power.t.test(n = floor(75 * (1-dropout)), delta = 1.92, sd = 3.84)
