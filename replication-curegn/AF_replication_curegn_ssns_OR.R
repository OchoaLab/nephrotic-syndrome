library(tidyverse)
library(genio)
library(qqman)
library(popgeninfer)
dir = '/datacommons/ochoalab/ssns_gwas/replication/saige_results/'
setwd( dir )

#############################
disease_subtype = "ssns_ctrl"
ancestry = "all"

# step 1: extract gnomad SNP id's
gnomad_AC <- read.table(paste0(dir, "gnomad/outfile/gnomad-genome_", disease_subtype, "_clean.txt"), sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
curegn_AC <- read.table(paste0(dir, "curegn/AF/ssns/curegn_AC_ssns.txt"), header = TRUE) 
cases_vs_control = merge(gnomad_AC, curegn_AC %>% mutate(CHROM = paste0("chr", CHR)), by = c("CHROM", "POS", "REF", "ALT")) 

# LRT on curegn vs Gnomad (case vs control)
# define function input
# x1/n1: curegn; x2/n2: gnomad
x1 <- cbind(cases_vs_control$curegn_AC_eur, cases_vs_control$curegn_AC_afr)
n1 <- cbind(cases_vs_control$curegn_AN_eur, cases_vs_control$curegn_AN_afr)
x2 <- cbind(cases_vs_control$AC_nfe, cases_vs_control$AC_afr)
n2 <- cbind(cases_vs_control$AN_nfe, cases_vs_control$AN_afr)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals_af_or <- af_test_or( x1, n1, x2, n2 )

# log10_pvalues_x <- -log10(pvals_af)
# log10_pvalues_y <- -log10(pvals_af_or)
# data <- data.frame(log10_pvalues_x, log10_pvalues_y)
# 
# ggplot(data, aes(x = log10_pvalues_x, y = log10_pvalues_y)) +
#   geom_point(alpha = 0.5) +
#   labs(x = "-log10 p-values (AF)", y = "-log10 p-values (AF OR)",
#        title = "Replication: CureGN vs gnomad Controls") + 
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
#   theme_minimal()
# 
cases_vs_control_test = cbind(cases_vs_control, pvals_af_or) 

chrom = c(6,10,16,6)
pos = c(32652506,28810849,11077745,32689478)
sig_table = cases_vs_control_test %>% mutate(CHROM = str_remove(CHROM, "chr")) %>%  filter(CHROM %in% chrom & POS %in% pos)

# suggestive sig
chrom = c(8, 7, 8)
pos = c(107771874, 3014907, 142322618)
sugg_sig_table = cases_vs_control_test %>% mutate(CHROM = str_remove(CHROM, "chr"))  %>% filter(CHROM %in% chrom & POS %in% pos)
