library(tidyverse)
library(genio)
library(qqman)
dir = '/datacommons/ochoalab/ssns_gwas/replication/saige_results/'
setwd( dir )

#############################
disease_subtype = "ns_ctrl"
ancestry = "sas"
#############################

# step 1: extract ukbb SNP id's
ukbb_AC <- read.csv(paste0(dir, "ukbb/output_file/UKBB_ns_ctrl_sas.csv"))
bristol_AC <- read.table(paste0(dir, "bristol/AF/",  disease_subtype,  "/", ancestry,  "/bristol_AC_ns_allage_sas.txt"), header = TRUE) %>% 
  mutate(CHROM = str_remove(CHROM, "chr"))

cases_vs_control = merge(ukbb_AC, bristol_AC, by.x = c("Chromosome", "Position", "A1", "A2"), 
                         by.y = c("CHROM", "POS", "REF", "ALT")) %>% 
  select(Chromosome, Position, A1, A2,   
         ukbb_AC_sas, ukbb_AN_sas,
         bristol_AC_sas, bristol_AN_sas) %>% distinct()

# LRT on Bristol vs Gnomad (case vs control)
# define function input
# x1/n1: bistrol; x2/n2: gnomad
x1 <- cbind(cases_vs_control$bristol_AC_sas)
n1 <- cbind(cases_vs_control$bristol_AN_sas)
x2 <- cbind(cases_vs_control$ukbb_AC_sas)
n2 <- cbind(cases_vs_control$ukbb_AN_sas)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals_b <- af_test_or( x1, n1, x2, n2 )

output_plot = cbind(cases_vs_control, pvals_b)
# write.table(output_plot, paste0(dir, "bristol/AF/annotated_results/replication_pvals_", disease_subtype, "_", ancestry, ".txt"), 
#             row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

chrom = c(6)
pos = c(32604403)
sig_table = output_plot %>% filter(Chromosome %in% chrom & Position %in% pos)
