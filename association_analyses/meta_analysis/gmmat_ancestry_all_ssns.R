library(GMMAT)
library(genio)
library(tidyverse)


name <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/ssns_ctr_mac20_ancestry.grm"
grm = as.matrix(read_grm(name)$kinship)

phenofile = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/pheno_ssns_ancestry.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) %>% mutate(V3 = ifelse(V3 == 1, 0, 1))
covarfile = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/covar_ssns_ancestry.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE)

eigenfile = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/ssns_ctr_mac20_ancestry.eigenvec", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 

pheno = merge(phenofile %>% dplyr::rename(trait = V3), covar_sub %>% dplyr::rename(sex = V3, race = V4), by = c("V1","V2"), all.y = TRUE)
pheno_all = merge(pheno, eigenfile, by = c("V1","V2")) %>% dplyr::rename(famid = V1, id = V2)

model0 <- glmmkin(trait ~  sex + race + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, data = pheno_all, kins = grm, id = "id",family = binomial(link = "logit"))


glmm.score(model0, infile = "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/ssns_ctr_mac20_ancestry", outfile ="glmm.score.bed.meta_ancestry_combined_ssns_ctr.txt")


