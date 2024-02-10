library(tidyverse)
library(popgeninfer)
# chr = c("chr6","chr6","chr6","chr6","chr6","chr10","chr6","chr17","chr5","chr22","chr11","chr17","chr16","chr18","chr13","chr8","chr17")
# pos = c(32667852,32684117,32620075,32542235,32658788,28810849,32589312,56125448,23079630,43217545,110737599,51793737,11077745,2727044,78236524,109898199, 47740256)
# ancestry = c("All", "All", "European", "African", "South Asian", "All", "African", "All", "African", "European", "All", "European", "All", "All", "All", "All", "All")
# ukbb_list = cbind(chr, pos, ancestry) %>% as.data.frame()

ukbb_controls = read.csv('/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ukbb/ssns_ukbb_controls.csv')[1:17, 1:13]
ukbb_all = ukbb_controls %>% filter(ancestry == "All" & !rsid %in% c("rs151109286", 'rs11213573')) %>% drop_na()

bristol_ssns = read.table('/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/bristol_AC_ssns_allage_allsnps.txt', header = TRUE)
#7075
# extract from bristol
control_snp = ukbb_all %>% separate(allele, c("A1", "A2")) %>% mutate(SNP = paste0(CHR, ":", POS, ":",  A2, ":", A1)) %>% select(SNP)
write.table(control_snp,
            "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ukbb/control_snp_id_ssns_all.txt", 
            row.names=FALSE, quote=FALSE, col.names = FALSE)

# read AF files
sas_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ukbb/bristolAC_sas.acount")
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
sas_ac_clean = sas_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristolAC_sas = ALT_FREQS, bristolAN_sas = OBS_CT)

eur_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ukbb/bristolAC_eur.acount")
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristolAC_eur = ALT_FREQS, bristolAN_eur = OBS_CT)

afr_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ukbb/bristolAC_afr.acount")
colnames(afr_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
afr_ac_clean = afr_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristolAC_afr = ALT_FREQS, bristolAN_afr = OBS_CT)

merge1 = merge(sas_ac_clean, eur_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)
merge2 = merge(merge1, afr_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)

case_vs_control_all = merge(merge2, ukbb_all, by.x = c("CHROM", "POS"), by.y = c("CHR", "POS")) %>% 
  mutate(UKBB_eur_AC = as.numeric(UKBB_eur_AC), UKBB_afr_AC = as.numeric(UKBB_afr_AC), UKBB_sas_AC = as.numeric(UKBB_sas_AC),
         UKBB_eur_AN = as.numeric(UKBB_eur_AN), UKBB_afr_AN = as.numeric(UKBB_afr_AN), UKBB_sas_AN = as.numeric(UKBB_sas_AN))

# AF test
x1 <- cbind(case_vs_control_all$bristolAC_eur,  case_vs_control_all$bristolAC_sas, case_vs_control_all$bristolAC_afr) 
n1 <- cbind(case_vs_control_all$bristolAN_eur, case_vs_control_all$bristolAN_sas, case_vs_control_all$bristolAN_afr)
x2 <- cbind(case_vs_control_all$UKBB_eur_AC, case_vs_control_all$UKBB_afr_AC, case_vs_control_all$UKBB_sas_AC) 
n2 <- cbind(case_vs_control_all$UKBB_eur_AN, case_vs_control_all$UKBB_afr_AN, case_vs_control_all$UKBB_sas_AN) 
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals <- af_test( x1, n1, x2, n2 )$pval
case_vs_control_all_af = cbind(case_vs_control_all, pvals)

# european
ukbb_eur = ukbb_controls %>% filter(ancestry == "European") %>% drop_na() %>% 
  select(-UKBB_afr_AC, -UKBB_afr_AN, -UKBB_sas_AC, -UKBB_sas_AN)
bristol_ssns = read.table('/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/eur/bristol_AC_ssns_eur.txt', header = TRUE)
case_vs_control_eur = merge(bristol_ssns, ukbb_eur, by.x = c("CHROM", "POS"), by.y = c("CHR", "POS")) %>% 
  select(CHROM, POS, REF, ALT, rsid, allele, type, gene, ancestry, everything())

# AF test
x1 <- case_vs_control_eur$bristol_AC_eur %>% as.numeric()
n1 <- case_vs_control_eur$bristol_AN_eur %>% as.numeric()
x2 <- case_vs_control_eur$UKBB_eur_AC %>% as.numeric()
n2 <- case_vs_control_eur$UKBB_eur_AN %>% as.numeric()
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals <- af_test( x1, n1, x2, n2 )$pval
case_vs_control_eur_af = cbind(case_vs_control_eur, pvals)
# african
ukbb_afr = ukbb_controls %>% filter(ancestry == "African") %>% drop_na() %>% 
  select(-UKBB_eur_AC, -UKBB_eur_AN, -UKBB_sas_AC, -UKBB_sas_AN)
# extract from bristol
control_snp = ukbb_afr %>% separate(allele, c("A1", "A2")) %>% mutate(SNP = paste0(CHR, ":", POS, ":",  A2, ":", A1)) %>% select(SNP)
write.table(control_snp,
            "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ukbb/control_snp_id_ssns_afr.txt", 
            row.names=FALSE, quote=FALSE, col.names = FALSE)
## low MAC in AFR

# south asian
ukbb_sas = ukbb_controls %>% filter(ancestry == "South Asian") %>% drop_na() %>% 
  select(-UKBB_eur_AC, -UKBB_eur_AN, -UKBB_afr_AC, -UKBB_afr_AN)
bristol_ssns = read.table('/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/sas/bristol_AC_ssns_sas.txt', header = TRUE)
# extract from bristol
control_snp = ukbb_sas %>% separate(allele, c("A1", "A2")) %>% mutate(SNP = paste0(CHR, ":", POS, ":",  A2, ":", A1)) %>% select(SNP)
write.table(control_snp,
            "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ukbb/control_snp_id_ssns_sas.txt", 
            row.names=FALSE, quote=FALSE, col.names = FALSE)

sas_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ukbb/bristolAC_sas_nofilter.acount")
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
sas_ac_clean = sas_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% 
  dplyr::rename(bristol_AC_sas = ALT_FREQS, bristol_AN_sas = OBS_CT) %>% 
  mutate(ALT = ifelse(ALT == "TRUE", "T", ALT))
case_vs_control_sas = merge(sas_ac_clean, ukbb_sas, by.x = c("CHROM", "POS"), by.y = c("CHR", "POS")) %>% 
  select(CHROM, POS, REF, ALT, rsid, allele, type, gene, ancestry, everything())

# AF test
x1 <- case_vs_control_sas$bristol_AC_sas %>% as.numeric()
n1 <- case_vs_control_sas$bristol_AN_sas %>% as.numeric()
x2 <- case_vs_control_sas$UKBB_sas_AC %>% as.numeric()
n2 <- case_vs_control_sas$UKBB_sas_AN %>% as.numeric()
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals <- af_test( x1, n1, x2, n2 )$pval
case_vs_control_sas_af = cbind(case_vs_control_sas, pvals)
