library(tidyverse)
library(genio)
library(qqman)
dir = '/datacommons/ochoalab/ssns_gwas/replication/saige_results/'
setwd( dir )

#############################
disease_subtype = "ssns_ctrl"
ancestry = "all"
#############################

# step 1: extract ukbb SNP id's
ukbb_AC <- read.csv(paste0(dir, "ukbb/output_file/UKBB_ssns_ctrl.csv"))
bristol_AC <- read.table(paste0(dir, "bristol/AF/",  disease_subtype, "/bristol_AC_ssns_allage.txt"), header = TRUE) %>% 
  mutate(CHROM = str_remove(CHROM, "chr"))

cases_vs_control = merge(ukbb_AC, bristol_AC, by.x = c("Chromosome", "Position", "A1", "A2"), 
                         by.y = c("CHROM", "POS", "REF", "ALT")) %>% 
  select(Chromosome, Position, A1, A2,  ukbb_AC_afr, ukbb_AN_afr, 
         ukbb_AC_eur, ukbb_AN_eur, ukbb_AC_sas, ukbb_AN_sas, bristol_AC_afr, bristol_AN_afr, 
         bristol_AC_eur, bristol_AN_eur, bristol_AC_sas, bristol_AN_sas) %>% distinct()

# LRT on Bristol vs Gnomad (case vs control)
# define function input
# x1/n1: bistrol; x2/n2: gnomad
x1 <- cbind(cases_vs_control$bristol_AC_sas, cases_vs_control$bristol_AC_eur, cases_vs_control$bristol_AC_afr)
n1 <- cbind(cases_vs_control$bristol_AN_sas, cases_vs_control$bristol_AN_eur, cases_vs_control$bristol_AN_afr)
x2 <- cbind(cases_vs_control$ukbb_AC_sas, cases_vs_control$ukbb_AC_eur, cases_vs_control$ukbb_AC_afr)
n2 <- cbind(cases_vs_control$ukbb_AN_sas, cases_vs_control$ukbb_AN_eur, cases_vs_control$ukbb_AN_afr)
# p-values for forward alignment (only alignment for non-revcomp SNPs)

pvals_af_or <- af_test_or( x1, n1, x2, n2 )

# log10_pvalues_x <- -log10(pvals_af)
# log10_pvalues_y <- -log10(pvals_af_or)
# data <- data.frame(log10_pvalues_x, log10_pvalues_y)
# 
# ggplot(data, aes(x = log10_pvalues_x, y = log10_pvalues_y)) +
#   geom_point(alpha = 0.5) +
#   labs(x = "-log10 p-values (AF)", y = "-log10 p-values (AF OR)",
#        title = "Replication: Bristol vs UKBB Controls") + 
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
#   theme_minimal()

# hist(pvals_b, breaks = 50, main = "Bristol vs UKBB")
# qq(pvals_b, main = "Bristol vs UKBB")


output_plot = cbind(cases_vs_control, pvals_af_or)
# write.table(output_plot, paste0(dir, "bristol/AF/annotated_results/replication_pvals_", disease_subtype, ".txt"), 
#             row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

chrom = c(16, 6, 6, 6)
pos = c(11077745, 32689478, 32669238, 32658788)
sig_table = output_plot %>% filter(Chromosome %in% chrom & Position %in% pos)

# suggestive sig
chrom = c(17, 8, 7, 8)
pos = c(56125448, 107771874, 3014907, 142322618)
sugg_sig_table = output_plot %>% filter(Chromosome %in% chrom & Position %in% pos)
