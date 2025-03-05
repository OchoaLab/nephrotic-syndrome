library(tidyverse)
library(genio)
library(qqman)
library(popgeninfer)
dir = '/datacommons/ochoalab/ssns_gwas/replication/saige_results/'
setwd( dir )

#############################
disease_subtype = "srns_ctrl"
ancestry = "all"


# step 1: extract gnomad SNP id's
gnomad_AC <- read.table(paste0(dir, "gnomad/outfile/gnomad-genome_", disease_subtype, "_clean.txt"), sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
curegn_AC <- read.table(paste0(dir, "curegn/AF/srns/curegn_AC_srns.txt"), header = TRUE) 
cases_vs_control = merge(gnomad_AC, curegn_AC %>% mutate(CHROM = paste0("chr", CHR)), by = c("CHROM", "POS", "REF", "ALT")) 

# LRT on curegn vs Gnomad (case vs control)
# define function input
# x1/n1: bistrol; x2/n2: gnomad
x1 <- cbind(cases_vs_control$curegn_AC_eur, cases_vs_control$curegn_AC_afr)
n1 <- cbind(cases_vs_control$curegn_AN_eur, cases_vs_control$curegn_AN_afr)
x2 <- cbind(cases_vs_control$AC_nfe, cases_vs_control$AC_afr)
n2 <- cbind(cases_vs_control$AN_nfe, cases_vs_control$AN_afr)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals <- af_test_or( x1, n1, x2, n2 )

cases_vs_control_test = cbind(cases_vs_control, pvals) 

chrom = c(22,22,22,22,22,22)
pos = c(36267167, 36260509,36260398,36271108,36265860,36265988)
sig_table = cases_vs_control_test %>% mutate(CHROM = str_remove(CHROM, "chr")) %>%  filter(CHROM %in% chrom & POS %in% pos)
