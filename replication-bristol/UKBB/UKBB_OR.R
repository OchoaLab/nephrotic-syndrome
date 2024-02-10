library(tidyverse)
library(epitools)
#### OR calculations
# manuscript table
chr = c("chr6","chr6","chr6","chr6","chr6","chr10","chr6","chr17","chr5","chr22","chr11","chr17","chr16","chr18","chr13","chr8","chr17")
pos = c(32667852,32684117,32620075,32542235,32658788,28810849,32589312,56125448,23079630,43217545,110737599,51793737,11077745,2727044,78236524,109898199, 47740256)
ancestry = c("All", "All", "European", "African", "South Asian", "All", "African", "All", "African", "European", "All", "European", "All", "All", "All", "All", "All")
ukbb_list = cbind(chr, pos, ancestry) %>% as.data.frame()
chr_ = str_remove(chr, "chr")
# discovery analysis snps
## european
ukbb_eur = ukbb_list %>% filter(ancestry == "European") %>% mutate(chr = str_remove(chr, "chr"))
ssns_eur_results = read.table('/datacommons/ochoalab/ssns_gwas/imputed/annotate/locuszoom/ssns_ctrl_eur.txt', header = TRUE)  %>% 
  filter(CHR %in% ukbb_eur$chr & POS %in% ukbb_eur$pos)

write.table(ssns_eur_results %>% pull(SNP), paste0("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/eur/snp_id_from_table.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE) 

## sas
ukbb_sas = ukbb_list %>% filter(ancestry == "South Asian") %>% mutate(chr = str_remove(chr, "chr"))
ssns_sas_results = read.table('/datacommons/ochoalab/ssns_gwas/imputed/annotate/locuszoom/ssns_ctrl_sas.txt', header = TRUE)  %>% 
  filter(CHR %in% ukbb_sas$chr & POS %in% ukbb_sas$pos)
write.table(ssns_sas_results %>% pull(SNP), paste0("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/sas/snp_id_from_table.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE) 


## afr
ukbb_afr = ukbb_list %>% filter(ancestry == "African") %>% mutate(chr = str_remove(chr, "chr"))
ssns_afr_results = read.table('/datacommons/ochoalab/ssns_gwas/imputed/annotate/locuszoom/ssns_ctrl_afr.txt', header = TRUE)  %>% 
  filter(CHR %in% ukbb_afr$chr & POS %in% ukbb_afr$pos)

write.table(ssns_afr_results %>% pull(SNP), paste0("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/afr/snp_id_from_table.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE) 

## all ancestry
ukbb_all = ukbb_list %>% filter(ancestry == "All") %>% mutate(chr = str_remove(chr, "chr"))
ssns_all_results = read.table('/datacommons/ochoalab/ssns_gwas/imputed/annotate/locuszoom/ssns_ctrl_all.txt', header = TRUE)  %>% 
  filter(CHR %in% ukbb_all$chr & POS %in% ukbb_all$pos)

write.table(ssns_all_results %>% pull(SNP), paste0("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/snp_id_from_table.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE) 

# extract AF from Bristol
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --extract $dir/ssns_ctrl/snp_id_from_table.txt --freq 'counts' --out $dir/ssns_ctrl/bristolAC_ukbb
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_sas_allage.txt --extract $dir/ssns_ctrl/sas/snp_id_from_table.txt --freq 'counts' --out $dir/ssns_ctrl/sas/bristolAC_ukbb
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_af_allage.txt --extract $dir/ssns_ctrl/afr/snp_id_from_table.txt --freq 'counts' --out $dir/ssns_ctrl/afr/bristolAC_ukbb
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_eur_allage.txt --extract $dir/ssns_ctrl/eur/snp_id_from_table.txt --freq 'counts' --out $dir/ssns_ctrl/eur/bristolAC_ukbb

# 
eur_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/eur/bristolAC_ukbb.acount")
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "Bristol_AC", "Bristol_N")
sas_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/sas/bristolAC_ukbb.acount")
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "Bristol_AC", "Bristol_N")
afr_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/afr/bristolAC_ukbb.acount")
colnames(afr_ac) = c("CHROM", "ID", "REF", "ALT", "Bristol_AC", "Bristol_N")
all_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/bristolAC_ukbb.acount")
colnames(all_ac) = c("CHROM", "ID", "REF", "ALT", "Bristol_AC", "Bristol_N")


# eur
UKBB_AC = c(545872, 37970, 660483)
UKBB_N = c(899544, 917616, 916894)
eur = cbind(eur_ac, UKBB_AC, UKBB_N)

# sas
UKBB_AC = 4452
UKBB_N = 19336
sas = cbind(sas_ac, UKBB_AC, UKBB_N)

# afr
UKBB_AC = 271
UKBB_N = 17526
afr = cbind(afr_ac, UKBB_AC, UKBB_N)

# all
UKBB_AC = c(27962, NA, 322544, 23148)
UKBB_N = c(930608, NA, 980616, 981074)
all = cbind(all_ac, UKBB_AC, UKBB_N)

AC_test = gdata::combine(eur, sas, afr, all) %>% drop_na()
# variables to be feed into oddsratio.wald() function
a = AC_test$Bristol_AC
b = AC_test$Bristol_N - AC_test$Bristol_AC
c = AC_test$UKBB_AC
d = AC_test$UKBB_N - AC_test$UKBB_AC

lower_CI = exp(log(Inf) - 1.96*sqrt(1/a + 1/b + 1/c + 1/d))
upper_CI = exp(log(Inf) + 1.96*sqrt(1/a + 1/b + 1/c + 1/d))

out = data.frame()
fisher_exact = c()
for (i in 1:length(a)){
  # create a freq_table for each entry for num_risk_allele
  freq_table = matrix(c(a[i], b[i], c[i], d[i]), ncol = 2)
  # if there is NA in the freq_table, ignore (preprocessing should not allow this to happen)
  if (any(is.na(freq_table)) == TRUE){
    next
  } else {
    OR_out = oddsratio.wald(freq_table)$measure %>% data.frame() %>%
      rename(OR = estimate, lower_CI = lower, upper_CI = upper) %>%
      mutate(number_risk_allele = as.factor(paste(i-1))) %>% select(number_risk_allele, everything())
    rownames(OR_out) <- NULL
    fisher = oddsratio.wald(freq_table)$p.value %>% as.data.frame() %>% pull(fisher.exact)
    fisher_exact = c(fisher_exact, fisher[2])
    out = rbind(out, OR_out[2,])
  }
  output = cbind(out, fisher_exact) %>%
    rename(num_risk_allele = number_risk_allele)
}



