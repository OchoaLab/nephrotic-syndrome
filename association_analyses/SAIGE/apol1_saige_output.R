library(tidyverse)
output_ssns = read.table('/datacommons/ochoalab/ssns_gwas/saige/APOL1/ssns_apol1.txt', header = TRUE)
output_srns = read.table('/datacommons/ochoalab/ssns_gwas/saige/APOL1/srns_apol1.txt', header = TRUE)

# calculate OR, CI
calculate_OR_and_CI <- function(beta, se) {
  OR <- exp(beta)
  lower_bound <- exp(beta - 1.96 * se)
  upper_bound <- exp(beta + 1.96 * se)
  return(list(OR = OR, lower_bound = lower_bound, upper_bound = upper_bound))
}

calculate_OR_and_CI(output_ssns$BETA, output_ssns$SE) # 1.04(0.55-1.65)
output_ssns$p.value # 0.862
output_ssns$BETA # 0.0405
output_ssns$SE # 0.233

calculate_OR_and_CI(output_srns$BETA, output_srns$SE) # 2.62 (1.31, 5.25)

output_srns$p.value # 0.006629751
output_srns$BETA # 0.96
output_srns$SE # 0.355

