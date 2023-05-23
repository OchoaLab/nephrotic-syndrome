library(tidyverse)

### ssns ###
ssns_ctr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_cond1 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_cond1_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_cond2 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_cond2_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_cond3 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_cond3.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_cond3 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_cond3_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
            
ssns_ctr_cond4 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_cond4.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_cond4 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_cond4_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
            
ssns_ctr_afr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_afr %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_afr_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_afr_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_afr_cond1 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_afr_cond1_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_afr_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_afr_cond2 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_afr_cond2_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_eur <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_eur %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_eur_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_eur_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_eur_cond1 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_eur_cond1_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_eur_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_eur_cond2 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_eur_cond2_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_sas <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_sas %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_sas_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_sas_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_sas_cond1 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_sas_cond1_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ssns_ctr_sas_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ssns_ctr_sas_cond2 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ssns_ctr_sas_cond2_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

### ns ###
ns_ctr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)


ns_ctr_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_cond1 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_cond1_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)


ns_ctr_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_cond2 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_cond2_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ns_ctr_cond3 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_cond3.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_cond3 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_cond3_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ns_ctr_afr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_afr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_afr %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_afr_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)


ns_ctr_afr_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_afr_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_afr_cond1 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_afr_cond1_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)


ns_ctr_eur <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_euro_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_eur %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_eur_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)


ns_ctr_eur_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_euro_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_eur_cond1 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_eur_cond1_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)


ns_ctr_sas <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_sas_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_sas %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_sas_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

ns_ctr_sas_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_sas_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(ns_ctr_sas_cond1 %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/ns_ctr_sas_cond1_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

### srns ###

srns_ctr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/srns_ctr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(srns_ctr %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/srns_ctr_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

srns_ssns <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/srns_ssns_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(srns_ssns %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/srns_ssns_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

### meta-analysis ###
meta_ssns <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_meta.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(meta_ssns %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/meta_ssns_ctr_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

meta_ns <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_meta.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
write.table(meta_ns %>% select(SNP),
            "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/snp_id/meta_ns_ctr_snpid.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)