library(tidyverse)

sas_ac <- read.table("b_ssns_sas_allage.acount") 
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
sas_ac_clean = sas_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% rename(bristol_AC_sas = ALT_FREQS, bristol_AN_sas = OBS_CT)

euro_ac <- read.table("b_ssns_euro_allage.acount") 
colnames(euro_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
euro_ac_clean = euro_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% rename(bristol_AC_euro = ALT_FREQS, bristol_AN_euro = OBS_CT)

af_ac <- read.table("b_ssns_af_allage.acount") 
colnames(af_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
af_ac_clean = af_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% rename(bristol_AC_af = ALT_FREQS, bristol_AN_af = OBS_CT)

merge1 = merge(sas_ac_clean, euro_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)
merge2 = merge(merge1, af_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)

write.table(merge2, 
            "bristol_AC_ssns_allage.txt", 
            row.names=FALSE, quote=FALSE, col.names = TRUE)
