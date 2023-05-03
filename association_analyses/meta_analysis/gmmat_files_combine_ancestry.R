
imp_b = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/imp_b.txt")
imp_sa = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/imp_sa.txt")
imp_w = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/imp_w.txt")

imp_all = rbind(imp_b, imp_sa, imp_w)

write.table(imp_all, "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/imp_all_ancestry.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)


pheno_b = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/b_pheno.txt")
pheno_sa = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/sa_pheno.txt")
pheno_w = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/w_pheno.txt")

pheno_all = rbind(pheno_b, pheno_sa, pheno_w)

write.table(pheno_all, "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/pheno_all_ancestry.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

covar_b = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/maf/b_covar.txt")
covar_sa = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/maf/sa_covar.txt")
covar_w = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/maf/w_covar.txt")

covar_all = rbind(covar_b, covar_sa, covar_w)

write.table(covar_all, "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/covar_all_ancestry.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)


## SSNS
# pheno
pheno_b = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_b.txt")
pheno_sa = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_sa.txt")
pheno_w = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_w.txt")

pheno_all = rbind(pheno_b, pheno_sa, pheno_w)

write.table(pheno_all, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/pheno_ssns_ancestry.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

covar_b = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_covar_b.txt")
covar_sa = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_covar_sa.txt")
covar_w = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_covar_w.txt")

covar_all = rbind(covar_b, covar_sa, covar_w)

write.table(covar_all, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/covar_ssns_ancestry.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(covar_all %>% select(V1, V2), "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/snpid_ssns_ancestry.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

# ssns_b = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/imp_b.txt")
# ssns_sa = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/imp_sa.txt")
# ssns_w = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/imp_w.txt")
# 
