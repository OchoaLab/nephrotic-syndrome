library(tidyverse)

sas_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_new/ssns_ctr_sas_gnomadsnps.acount") 
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
sas_ac_clean = sas_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% rename(disc_AC_sas = ALT_FREQS, disc_AN_sas = OBS_CT)

euro_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_new/ssns_ctr_euro_gnomadsnps.acount") 
colnames(euro_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
euro_ac_clean = euro_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% rename(disc_AC_euro = ALT_FREQS, disc_AN_euro = OBS_CT)

af_ac <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_new/ssns_ctr_afr_gnomadsnps.acount") 
colnames(af_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
af_ac_clean = af_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% rename(disc_AC_af = ALT_FREQS, disc_AN_af = OBS_CT)

merge1 = merge(sas_ac_clean, euro_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)
merge2 = merge(merge1, af_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)

write.table(merge2, 
            "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_new/discovery_AC_ssns_allage.txt", 
            row.names=FALSE, quote=FALSE, col.names = TRUE)
