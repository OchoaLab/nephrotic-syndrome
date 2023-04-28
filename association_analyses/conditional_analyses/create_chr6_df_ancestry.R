library(genio)
library(tidyverse)
library(BEDMatrix)

afr_mac <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/glmm.score.bed.ancestry_ssns_ctr_afr_mac20_PC_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
afr_mac_chr6 = afr_mac %>% filter(CHR == 6) %>% arrange(PVAL) %>% slice(1)
print(afr_mac_chr6$SNP)

# african ancestry
name <- '/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/ssns_ctr_afr_mac20'
ssns_bed = BEDMatrix(name, simple_names = TRUE)
chr6_geno = ssns_bed[,print(afr_mac_chr6$SNP)]
covar = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_covar_b.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE)
covar_chr6 = cbind(covar, chr6_geno)
write.table(covar_chr6, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/afr/covar_chr6_ssns_control_afr.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)



euro_mac <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/euro/glmm.score.bed.ancestry_ssns_ctr_euro_PC_mac20_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
euro_mac_chr6 = euro_mac %>% filter(CHR == 6) %>% arrange(PVAL) %>% slice(1)
print(euro_mac_chr6$SNP)
# european ancestry
name <- '/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/euro/ssns_ctr_euro_mac20'
ssns_bed = BEDMatrix(name, simple_names = TRUE)
chr6_geno = ssns_bed[,print(euro_mac_chr6$SNP)]
covar = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_covar_w.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE)
covar_chr6 = cbind(covar, chr6_geno)
write.table(covar_chr6, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/euro/covar_chr6_ssns_control_euro.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)



sa_mac <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/sas/glmm.score.bed.ancestry_ssns_ctr_sas_mac20_PC_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
sa_mac_chr6 = sa_mac %>% filter(CHR == 6) %>% arrange(PVAL) %>% slice(1)
print(sa_mac_chr6$SNP)
# south asian ancestry
name <- '/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/sas/ssns_ctr_sas_mac20'
ssns_bed = BEDMatrix(name, simple_names = TRUE)
chr6_geno = ssns_bed[,print(sa_mac_chr6$SNP)]
covar = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_covar_sa.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE)
covar_chr6 = cbind(covar, chr6_geno)
write.table(covar_chr6, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/sas/covar_chr6_ssns_control_sas.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)