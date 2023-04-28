library(genio)
library(tidyverse)
library(BEDMatrix)

gmmat_mac <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/glmm.score.bed.all_ssns_ctr_mac20_PC_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

gmmat_mac_chr6 = gmmat_mac %>% filter(CHR == 6) %>% arrange(PVAL) %>% slice(1)
print(gmmat_mac_chr6$SNP)

# full data: ssns vs control
name <- '/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20'
ssns_bed = BEDMatrix(name, simple_names = TRUE)
chr6_geno = ssns_bed[,print(gmmat_mac_chr6$SNP)]
covar = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_covar.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE)
covar_chr6 = cbind(covar, chr6_geno)
write.table(covar_chr6, "covar_chr6_ssns_control.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)



