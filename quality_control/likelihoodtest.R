##Likelihood ratio test, results analyzed in hwe_analysis.Rmd
library(tidyverse)

afreq_ssns_b <- read.table("/datacommons/ochoalab/tiffany_data/DATA_NS_MERGE/datasubset/ssns_control_black.afreq", sep = "\t")
colnames(afreq_ssns_b) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")

afreq_ssns_w <- read.table("/datacommons/ochoalab/tiffany_data/DATA_NS_MERGE/datasubset/ssns_control_white.afreq", sep = "\t")
colnames(afreq_ssns_w) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")

afreq_ssns_a <- read.table("/datacommons/ochoalab/tiffany_data/DATA_NS_MERGE/datasubset/ssns_control_asian.afreq", sep = "\t")
colnames(afreq_ssns_a) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")

tgp_asian <- read.table("/datacommons/ochoalab/tiffany_data/DATA_NS_MERGE/datasubset/tgp_asian.afreq", sep = "\t")
colnames(tgp_asian) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")

tgp_white <- read.table("/datacommons/ochoalab/tiffany_data/DATA_NS_MERGE/datasubset/tgp_white.afreq", sep = "\t")
colnames(tgp_white) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")

tgp_black <- read.table("/datacommons/ochoalab/tiffany_data/DATA_NS_MERGE/datasubset/tgp_black.afreq", sep = "\t")
colnames(tgp_black) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")

#black
## define variables
ssns_b = afreq_ssns_b %>% mutate(alt_allele_count = ALT_FREQS * OBS_CT)
x1 = ssns_b$alt_allele_count
n1 = ssns_b$OBS_CT

tgp_b = tgp_black %>% mutate(alt_allele_count = ALT_FREQS * OBS_CT)
x2 = tgp_b$alt_allele_count
n2 = tgp_b$OBS_CT

## alt pmf
p1 = x1/n1
pmf1 = dbinom(x1, n1, p1)
p2 = x2/n2
pmf2 = dbinom(x2, n2, p2)
altpmf = pmf1*pmf2

## null pmf
nullpmf = (x1+x2)/(n1+n2)

## lambda - converges asymptotically to chi-sqaured distribution
lamdba_lr = -2*log(nullpmf/altpmf)
output = pchisq(q = lamdba_lr, df = 1, lower.tail = FALSE)

saveRDS(output, file = "chisq_lr_black.rds")

#asian
## define variables
ssns_a = afreq_ssns_a %>% mutate(alt_allele_count = ALT_FREQS * OBS_CT)
x1 = ssns_a$alt_allele_count
n1 = ssns_a$OBS_CT

tgp_a = tgp_asian %>% mutate(alt_allele_count = ALT_FREQS * OBS_CT)
x2 = tgp_a$alt_allele_count
n2 = tgp_a$OBS_CT

## alt pmf
p1 = x1/n1
pmf1 = dbinom(x1, n1, p1)
p2 = x2/n2
pmf2 = dbinom(x2, n2, p2)
altpmf = pmf1*pmf2

## null pmf
nullpmf = (x1+x2)/(n1+n2)

## lambda - converges asymptotically to chi-sqaured distribution
lamdba_lr = -2*log(nullpmf/altpmf)
output = pchisq(q = lamdba_lr, df = 1, lower.tail = FALSE)

saveRDS(output, file = "chisq_lr_asian.rds")

#white
## define variables
ssns_w = afreq_ssns_w %>% mutate(alt_allele_count = ALT_FREQS * OBS_CT)
x1 = ssns_w$alt_allele_count
n1 = ssns_w$OBS_CT

tgp_w = tgp_white%>% mutate(alt_allele_count = ALT_FREQS * OBS_CT)
x2 = tgp_w$alt_allele_count
n2 = tgp_w$OBS_CT

## alt pmf
p1 = x1/n1
pmf1 = dbinom(x1, n1, p1)
p2 = x2/n2
pmf2 = dbinom(x2, n2, p2)
altpmf = pmf1*pmf2

## null pmf
nullpmf = (x1+x2)/(n1+n2)

## lambda - converges asymptotically to chi-sqaured distribution
lamdba_lr = -2*log(nullpmf/altpmf)
output = pchisq(q = lamdba_lr, df = 1, lower.tail = FALSE)

saveRDS(output, file = "chisq_lr_white.rds")