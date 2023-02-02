library(GMMAT)
library(genio)
library(tidyverse)

# update covar file, identify top SNP
b_mac <- read.table("../glmm.score.bed.ancestry_black_ssns_ctr_mac20_PC.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

b_mac_chr6 = b_mac %>% filter(CHR == 6) %>% arrange(PVAL) %>% slice(1)
print(b_mac_chr6$SNP)

# african ancestry
name <- '../ssns_ctr_b_hwe_mac20'
plink = read_plink( name )
chr6_geno = plink$X[b_mac_chr6$SNP,]
covar = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_covar_b.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE)
covar_chr6 = cbind(covar, chr6_geno)
write.table(covar_chr6, "covar_chr6_b.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

# run GMMAT
name <- "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_b.grm"
grm = as.matrix(read_grm(name)$kinship)

phenofile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_b.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) %>% mutate(V3 = ifelse(V3 == 1, 0, 1))

covarfile = read.table("covar_chr6_b.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 

eigenfile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_b.eigenvec", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 

pheno = merge(phenofile %>% rename(trait = V3), covarfile %>% rename(sex = V3, chr6 = V4), by = c("V1","V2")) 
pheno_all = merge(pheno, eigenfile, by = c("V1","V2")) %>% rename(famid = V1, id = V2)

model0 <- glmmkin(trait ~  sex + chr6 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, data = pheno_all, kins = grm, id = "id",family = binomial(link = "logit"))
glmm.score(model0, infile = "../ssns_ctr_b_hwe_mac20", outfile ="glmm.score.bed.ancestry_black_ssns_ctr_mac20_PC_chr6.txt")


