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
x1 = ssns_b$alt_allele_count
n1 = ssns_b$OBS_CT
x2 = tgp_b$alt_allele_count
n2 = tgp_b$OBS_CT

likelihoodratio <- function(x1, n1, x2, n2, df){
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
  output_p = pchisq(q = lamdba_lr, df = df, lower.tail = FALSE)
  return(output_p)
}

output = likelihoodratio(x1, n1, x2, n2, 1)
# plot output: 
#hist(output, breaks = 50, main = "p-vals for White ancestry")
#qq(output, main = "QQ-plot for White ancestry")

## likelihood ratio test for 3 ancestry combined
likelihoodratio_3 = function(x1, n1, x2, n2, x3, n3, x4, n4, x5, n5, x6, n6){
  ## alt pmf
  p1 = x1/n1
  p2 = x2/n2
  pmf1 = dbinom(x1, n1, p1, log = TRUE)
  pmf2 = dbinom(x2, n2, p2, log = TRUE)
  
  p3 = x3/n3
  p4 = x4/n4
  pmf3 = dbinom(x3, n3, p3, log = TRUE)
  pmf4 = dbinom(x4, n4, p4, log = TRUE)
  
  p5 = x5/n5
  p6 = x6/n6
  pmf5 = dbinom(x5, n5, p5, log = TRUE)
  pmf6 = dbinom(x6, n6, p6, log = TRUE)
  
  altpmf_b = pmf1 + pmf2 
  altpmf_a = pmf3 + pmf4 
  altpmf_w = pmf5 + pmf6
  
  ## null pmf
  x_b = x1+x2
  x_a = x3+x4
  x_w = x5+x6
  n_b = n1+n2
  n_a = n3+n4
  n_w = n5+n6
  p_b = x_b/n_b
  p_a = x_a/n_a
  p_w = x_w/n_w
  
  nullpmf_b = dbinom(x1, n1, p_b, log = TRUE) + dbinom(x2, n2, p_b, log = TRUE)
  nullpmf_a = dbinom(x3, n3, p_a, log = TRUE) + dbinom(x4, n4, p_a, log = TRUE)
  nullpmf_w = dbinom(x5, n5, p_w, log = TRUE) + dbinom(x6, n6, p_w, log = TRUE)
  
  #likelihood ratio
  lamdba_lr = -2*((nullpmf_b - altpmf_b) + (nullpmf_a - altpmf_a) + (nullpmf_w - altpmf_w)) 
  # converge to chi-square distribution
  output_p = pchisq(q = lamdba_lr, df = 3, lower.tail = FALSE)
}
