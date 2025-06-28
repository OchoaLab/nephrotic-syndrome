library(tidyverse)
library(genio)
library(qqman)
library(popgeninfer)
dir = '/datacommons/ochoalab/ssns_gwas/replication/saige_results/'
setwd( dir )

#############################
disease_subtype = "ssns_ctrl"
ancestry = "eur"
#############################

# step 1: extract gnomad SNP id's
gnomad_AC <- read.table(paste0(dir, "gnomad/outfile/gnomad-genome_ssns_ctrl_eur_clean.txt"), sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
eur_ac <- read.table(paste0(dir, "curegn/AF/ssns/eur/curegnAC_eur_sub.acount")) 
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHROM", "POS"))  %>% dplyr::rename(curegn_AC_eur = ALT_FREQS, curegn_AN_eur = OBS_CT)

cases_vs_control = merge(gnomad_AC, eur_ac_clean %>% mutate(CHROM = paste0("chr", CHROM)), by = c("CHROM", "POS", "REF", "ALT"))

# LRT on curegn vs Gnomad (case vs control)
# define function input
# x1/n1: curegn; x2/n2: gnomad
x1 <- cbind(cases_vs_control$curegn_AC_eur)
n1 <- cbind(cases_vs_control$curegn_AN_eur)
x2 <- cbind(cases_vs_control$AC_nfe)
n2 <- cbind(cases_vs_control$AN_nfe)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals <- af_test_or( x1, n1, x2, n2 )
# define p-value threshold
cases_vs_control_test = cbind(cases_vs_control, pvals)

chrom = c(6)
pos = c(32669238)
sig_table = cases_vs_control_test %>% mutate(CHROM = str_remove(CHROM, "chr")) %>%  filter(CHROM %in% chrom & POS %in% pos)

