library(tidyverse)

# ### ssns vs control ### 
# afr_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/glmm.score.bed.ancestry_ssns_ctr_afr_mac20_PC_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# afr_original_sig = afr_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1) 
# write.table(afr_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
# euro_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/euro/glmm.score.bed.ancestry_ssns_ctr_euro_PC_mac20_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# euro_original_sig = euro_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL))# %>% select(CHR, POS, A2, A1)
# write.table(euro_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
# sas_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/sas/glmm.score.bed.ancestry_ssns_ctr_sas_mac20_PC_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# sas_original_sig = sas_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
# write.table(sas_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
# ssns_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/glmm.score.bed.all_ssns_ctr_mac20_PC_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# ssns_original_sig = ssns_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
# write.table(ssns_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
# ### ns vs control ### 
# afr_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/afr/glmm.score.bed.ancestry_afr_PC_mac20_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# afr_original_sig = afr_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
# write.table(afr_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_afr_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
# euro_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/euro/glmm.score.bed.ancestry_ns_euro_PC_mac20_no_hwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# euro_original_sig = euro_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
# write.table(euro_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_euro_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
# sas_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/sas/glmm.score.bed.ancestry_ns_sa_PC_mac20_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# sas_original_sig = sas_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
# write.table(sas_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_sas_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
# ns_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/glmm.score.bed.allPC_mac20_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# ns_original_sig = ns_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL))# %>% select(CHR, POS, A2, A1)
# write.table(ns_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
# ### srns vs control ###
 srns_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/glmm.score.bed.srns_ctr_match_mac20_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
 srns_original_sig = srns_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)\
 srns_original_sig %>% arrange(PVAL)
write.table(srns_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/srns_ctr_suggsig_summstat.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
# 
srns_ssns = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/srns_ssns/glmm.score.bed.srns_ssns_PC_mac20_nohwe.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
srns_ssns_sig = srns_ssns %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
srns_ssns_sig%>% arrange(PVAL)
write.table(srns_ssns_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/srns_ssns_suggsig_summstat_filter.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

## conditional analysis
### ssns vs control ### conditional 1
afr_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/afr/glmm.score.bed.ancestry_black_ssns_ctr_mac20_PC_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
afr_original_sig = afr_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1) 
write.table(afr_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat_cond1.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

euro_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/euro/glmm.score.bed.ancestry_white_ssns_ctr_20_PC_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
euro_original_sig = euro_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL))# %>% select(CHR, POS, A2, A1)
write.table(euro_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat_cond1.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

sas_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/sas/glmm.score.bed.ancestry_sa_ssns_ctr_20_PC_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
sas_original_sig = sas_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
write.table(sas_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat_cond1.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

ssns_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/glmm.score.bed.all_ssns_ctr_mac20_PC_nohwe_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_original_sig = ssns_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
write.table(ssns_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_cond1.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

### ssns vs control ### conditional 2
afr_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/conditional2/afr/glmm.score.bed.ancestry_black_ssns_ctr_mac20_PC_chr6_2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
afr_original_sig = afr_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1) 
write.table(afr_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat_cond2.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

euro_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/conditional2/euro/glmm.score.bed.ancestry_white_ssns_ctr_20_PC_chr6_2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
euro_original_sig = euro_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL))# %>% select(CHR, POS, A2, A1)
write.table(euro_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat_cond2.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

sas_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/conditional2/sas/glmm.score.bed.ancestry_sa_ssns_ctr_20_PC_chr6_2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
sas_original_sig = sas_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
write.table(sas_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat_cond2.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

ssns_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/conditional/conditional2/glmm.score.bed.all_ssns_ctr_mac20_PC_nohwe_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_original_sig = ssns_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
write.table(ssns_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_cond2.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

### ns vs control ### conditional 1
afr_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/conditional/afr/glmm.score.bed.ancestry_black_NS_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
afr_original_sig = afr_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
write.table(afr_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_afr_suggsig_summstat_cond1.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

euro_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/conditional/euro/glmm.score.bed.ancestry_white_NS_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
euro_original_sig = euro_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
write.table(euro_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_euro_suggsig_summstat_cond1.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

sas_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/conditional/sas/glmm.score.bed.ancestry_sa_NS_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
sas_original_sig = sas_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL)) #%>% select(CHR, POS, A2, A1)
write.table(sas_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_sas_suggsig_summstat_cond1.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

ns_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/conditional/glmm.score.bed.all.NSvsControl.PC_20_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ns_original_sig = ns_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL))# %>% select(CHR, POS, A2, A1)
write.table(ns_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_cond1.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

### ns vs control ### conditional 2
ns_original = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/conditional/conditional2/glmm.score.bed.all.NSvsControl.PC_20_chr6_2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ns_original_sig = ns_original %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL))# %>% select(CHR, POS, A2, A1)
write.table(ns_original_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_cond2.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

### meta-analysis
## ssns vs control:
ssns_meta = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/gmmat_ancestry_ssns_meta.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_meta_sig = ssns_meta %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL))# %>% select(CHR, POS, A2, A1)
write.table(ssns_meta_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_meta.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

## ns vs control:
ns_meta = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/meta_analysis_gmmat/gmmat_ancestry_ns_meta.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ns_meta_sig = ssns_meta %>% filter(PVAL < 1e-5) %>% mutate(PVAL = as.numeric(PVAL))# %>% select(CHR, POS, A2, A1)
write.table(ns_meta_sig, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_meta.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")
