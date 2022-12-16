library(GMMAT)
library(genio)
library(tidyverse)


name <- "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/gcta/ssns_tgp_impute_mac"
grm = as.matrix(read_grm(name)$kinship)

phenofile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/gcta/ssns_tgp_pheno_imp.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) %>% mutate(V3 = ifelse(V3 == 1, 0, 1))

covarfile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/gcta/covar_ns_tgp_new.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE)

pheno = merge(phenofile, covarfile, by = c("V1","V2"))
colnames(pheno) = c("famid", "id", "trait", "sex", "race")

model0 <- glmmkin(trait ~  sex + race, data = pheno, kins = grm, id = "id",family = binomial(link = "logit"))
glmm.score(model0, infile = "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ssns_tgp_imp_mac_10", outfile ="glmm.score.bed.all.txt")
