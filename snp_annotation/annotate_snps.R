library(tidyverse)

#### Load analysis suggestive sig dataframes
### ssns vs control ### 
ssns_ctr_afr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_ctr_afr_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_ctr_afr_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_afr_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

ssns_ctr_eur <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_ctr_eur_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_ctr_eur_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_euro_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

ssns_ctr_sas <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_ctr_sas_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_ctr_sas_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_sas_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

ssns_ctr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_ctr_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ssns_ctr_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

### ns vs control ### 
ns_ctr_afr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_afr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ns_ctr_afr_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_afr_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

ns_ctr_eur <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_euro_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ns_ctr_eur_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_euro_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

ns_ctr_sas <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_sas_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ns_ctr_sas_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_sas_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

ns_ctr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ns_ctr_cond1 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_cond1.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ns_ctr_cond2 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_cond2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

### srns vs control ### 
srns_ctr <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/srns_ctr_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

### srns vs ssns ### 
srns_ssns <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/srns_ssns_suggsig_summstat.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

### meta-analysis ###
meta_ssns <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig_summstat_meta.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
meta_ns <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ns_ctr_suggsig_summstat_meta.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)


#### combine dataframes, extract SNP id's
library(gdata)
combined_results = combine(ssns_ctr_afr, ssns_ctr_afr_cond1, ssns_ctr_afr_cond2, ssns_ctr_eur, ssns_ctr_eur_cond1, ssns_ctr_eur_cond2,
                           ssns_ctr_sas, ssns_ctr_sas_cond1, ssns_ctr_sas_cond2, ssns_ctr, ssns_ctr_cond1, ssns_ctr_cond2,
                           ns_ctr_afr, ns_ctr_afr_cond1, ns_ctr_eur, ns_ctr_eur_cond1, ns_ctr_sas, ns_ctr_sas_cond1, ns_ctr, 
                           ns_ctr_cond1, ns_ctr_cond2, srns_ctr, srns_ssns, meta_ssns, meta_ns) 
# 46106 SNPs, 11102 unique SNPs
unique_snpid = combined_results %>% select(SNP) %>% distinct()
write.table(unique_snpid, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/suggsig_uniqueID.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

#### run subset_sig.q
library(genio)
name <- '/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/combined_suggsig'
bim <- read_bim( name )
bim %>% group_by(id) %>% filter(n()>1)
bim = bim %>% mutate(id = paste0(chr, ":", pos))

data = read_plink(name)
X = data$X
X[1:5, 1:5]
rownames(X) <- bim$id
write_plink( "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/combined_suggsig_recode", X, bim, data$fam )

#### run_convert_vcf_rsid.q
rsid_vcf <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/rsid_combined_suggsig_recode.vcf.gz", sep = "\t", stringsAsFactors=FALSE, quote = "", header = FALSE)
rsid_vcf_sub = rsid_vcf[,1:5] %>% dplyr::rename(CHR = V1, POS = V2, rsid = V3, A2 = V4, A1 = V5)

merge_rsid = merge(combined_results, rsid_vcf_sub, by = c("CHR", "POS", "A1", "A2"), all.x = TRUE) %>% 
  select(CHR, POS, rsid, everything()) 

write.table(merge_rsid, "combined_list_pval_rsid.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")

#### for SNPNexus
nexus_list = combined_results %>% 
  mutate(type = "Chromosome", Strand = 1) %>% 
  select(type, Id = CHR, Position = POS, Allele1 = A1, Allele2 = A2, Strand) %>% distinct() %>% tail(1102)
write.table(nexus_list, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/nexus_list_unique_tail1102.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

