library(tidyverse)

# # simulated data
# ssns_x = rbinom(100, 2000, 0.5)
# tgp_x = rbinom(100, 2000, 0.5)
# 
# n1 = 2000
# n2 = 2000
# p1 = ssns_x/n1
# p2 = tgp_x/n2
# pmf1 = dbinom(ssns_x, n1, p1, log = TRUE)
# pmf2 = dbinom(tgp_x, n2, p2, log = TRUE)
# altpmf = pmf1+pmf2
# 
# x = ssns_x + tgp_x
# n = n1 + n2
# p = x/n
# nullpmf = dbinom(ssns_x, n1, p, log = TRUE)+dbinom(tgp_x, n2, p, log = TRUE) %>% as.vector()
# 
# lamdba_lr = -2*(nullpmf - altpmf)
# output = pchisq(q = lamdba_lr, df = 1, lower.tail = FALSE)
# hist(output)
# qq(output)

# Examples of input variables: 
# x1 = ssns_b$alt_allele_count
# n1 = ssns_b$OBS_CT
# x2 = tgp_b$alt_allele_count
# n2 = tgp_b$OBS_CT

likelihoodratio <- function(x1, n1, x2, n2){
  ## alt pmf
  p1 = x1/n1
  p2 = x2/n2
  pmf1 = dbinom(x1, n1, p1, log = TRUE)
  pmf2 = dbinom(x2, n2, p2, log = TRUE)
  altpmf = pmf1+pmf2
  
  ## null pmf
  x = x1+x2
  n = n1+n2
  p = x/n
  nullpmf = dbinom(x1, n1, p, log = TRUE) + dbinom(x2, n2, p, log = TRUE)
  #likelihood ratio
  lamdba_lr = -2*(nullpmf - altpmf)
  # converge to chi-square distribution
  output_p = pchisq(q = lamdba_lr, df = 1, lower.tail = FALSE)
  return(output_p)
}

# plot output: 
# hist(output, breaks = 50, main = "p-vals for White ancestry")
# qq(output, main = "QQ-plot for White ancestry")
