library(GMMAT)
library(genio)
library(tidyverse)


name <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/ssns_ctr_afr_mac20.grm"
grm = as.matrix(read_grm(name)$kinship)

phenofile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_b.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) %>% mutate(V3 = ifelse(V3 == 1, 0, 1))

covarfile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_covar_b.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 

eigenfile = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/ssns_ctr_afr_mac20.eigenvec", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 

pheno = merge(phenofile %>% rename(trait = V3), covarfile %>% rename(sex = V3), by = c("V1","V2")) 
pheno_all = merge(pheno, eigenfile, by = c("V1","V2")) %>% rename(famid = V1, id = V2)
table(pheno_all$trait)
write.table(pheno_all, "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_black_pheno_all.txt", row.names=FALSE, quote=FALSE, col.names=TRUE)

model0 <- glmmkin(trait ~  sex + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, data = pheno_all, kins = grm, id = "id",family = binomial(link = "logit"))
glmm.score(model0, infile = "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/ssns_ctr_afr_mac20", outfile ="glmm.score.bed.ancestry_ssns_ctr_afr_mac20_PC_nohwe.txt")

