library(tidyverse)
library(genio)
library(qqman)
dir = '/datacommons/ochoalab/ssns_gwas/replication/saige_results/'
setwd( dir )

#############################
disease_subtype = "ssns_ctrl"
ancestry = "all"
#############################

# step 1: extract ukbb SNP id's
ukbb_AC <- read.csv(paste0(dir, "ukbb/output_file/UKBB_ssns_ctrl.csv"))

ukbb_snp = ukbb_AC %>% mutate(SNP = paste0("chr", Chromosome, ":", Position, ":",  A1, ":", A2)) %>% select(SNP)
write.table(ukbb_snp, paste0(dir, "bristol/AF/", disease_subtype, "/ukbb_control_snp_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)

# # step 2: run allele_freq_discovery.q on NS data to remove false positive SNPs
# # load phenotype and fixed covariates file
# # write id files to calculate AC per ancestry
data <- read_tsv( '/datacommons/ochoalab/ssns_gwas/imputed/patient-data.txt.gz', col_types = 'cccccdciiii' )
table(data$srns_ctrl, data$ancestry)
# sas_id = data %>% filter(diagnosis == "Control" & ancestry == "sas") %>% select(id)
# afr_id = data %>% filter(diagnosis == "Control" & ancestry == "afr") %>% select(id)
# eur_id = data %>% filter(diagnosis == "Control" & ancestry == "eur") %>% select(id)
# all_id = data %>% filter(diagnosis == "Control") %>% select(id)
# 
# write.table(cbind(0, sas_id), paste0(dir, disease_subtype, "/ctrl_sas_id.txt"),
#             row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, afr_id), paste0(dir, disease_subtype, "/ctrl_afr_id.txt"),
#             row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, eur_id), paste0(dir, disease_subtype, "/ctrl_eur_id.txt"),
#             row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, all_id), paste0(dir, disease_subtype, "/ctrl_all_id.txt"),
#             row.names=FALSE, quote=FALSE, col.names = FALSE)
# 
# # step 3: run allele_freq_discovery.q 
sas_ac <- read.table(paste0(dir, "bristol/AF/ssns_ctrl/discoveryAC_sas_ukbbsnps.acount"))
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
sas_ac_clean = sas_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(disc_AC_sas = ALT_FREQS, disc_AN_sas = OBS_CT)

eur_ac <- read.table(paste0(dir, "bristol/AF/ssns_ctrl/discoveryAC_eur_ukbbsnps.acount"))
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(disc_AC_eur = ALT_FREQS, disc_AN_eur = OBS_CT)

afr_ac <- read.table(paste0(dir, "bristol/AF/ssns_ctrl/discoveryAC_afr_ukbbsnps.acount"))
colnames(afr_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
afr_ac_clean = afr_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(disc_AC_afr = ALT_FREQS, disc_AN_afr = OBS_CT)

merge1 = merge(sas_ac_clean, eur_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)
merge2 = merge(merge1, afr_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE) 

write.table(merge2, paste0(dir, "bristol/AF/", disease_subtype, "/discovery_AC_ancestry_combined.txt"),
            row.names=FALSE, quote=FALSE, col.names = TRUE)

# step 4: Merge discovery AC and UKBB AC
discovery_AC <- read.table(paste0(dir, "bristol/AF/", disease_subtype, "/discovery_AC_ancestry_combined.txt"), 
                           header = TRUE, sep = " ")

control_vs_control = merge(ukbb_AC %>% mutate(CHROM = paste0("chr", Chromosome)) %>% 
                             select(-Chromosome) %>% dplyr::rename(POS = Position), discovery_AC, 
                           by.x = c("CHROM", "POS", "A1", "A2"), by.y = c("CHROM", "POS", "REF", "ALT")) %>% 
  select(CHROM, POS, A1, A2, rsid, ukbb_AC_afr, ukbb_AN_afr, ukbb_AC_eur, ukbb_AN_eur, ukbb_AC_sas, ukbb_AN_sas,
         disc_AC_afr, disc_AN_afr, disc_AC_eur, disc_AN_eur, disc_AC_sas, disc_AN_sas) %>% distinct()


# step 5: AF test for discovery vs UKBB
library(popgeninfer)
# define function input
# x1/n1: discovery set; x2/n2: UKBB
x1 <- cbind(control_vs_control$disc_AC_sas, control_vs_control$disc_AC_eur, control_vs_control$disc_AC_afr)
n1 <- cbind(control_vs_control$disc_AN_sas, control_vs_control$disc_AN_eur, control_vs_control$disc_AN_afr )
x2 <- cbind(control_vs_control$ukbb_AC_sas, control_vs_control$ukbb_AC_eur, control_vs_control$ukbb_AC_afr)
n2 <- cbind(control_vs_control$ukbb_AN_sas, control_vs_control$ukbb_AN_eur, control_vs_control$ukbb_AN_afr)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals <- af_test( x1, n1, x2, n2 )$pval

# define p-value threshold
samplesize = length(pvals)
pcut = 0.05/samplesize

# sort by p-values
control_test = cbind(control_vs_control, pvals) %>% arrange(pvals) %>% 
  dplyr::rename(CHR = CHROM, BP = POS, P = pvals) %>% mutate(CHR = as.numeric(str_remove(CHR, "chr")), SNP = paste0(CHR, ":", BP)) %>% 
  arrange(CHR, BP) %>% mutate(snp_test = ifelse(P < pcut, "imputation bias", "replication"))
#table(control_test$CHR, control_test$sig)
table(control_test$snp_test)

list_to_remove = cbind(control_vs_control, pvals) %>% arrange(pvals) %>% filter(pvals < pcut) 
list_to_remove

hist(pvals, breaks = 50, main = "discovery control vs ukbb control (afr/eur/sas)")
qq(pvals, main = "discovery control vs ukbb control")

# step 6: Bristol AF test with UKBB with clean set of SNPs: 
## write SNPs into file 
## write per ancestry id's
# # write subset id's with phenotype
pheno_ <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/admixture/pca1_pheno_df.txt", sep = "\t", header = TRUE) %>%
  select(IID, Sex, RACE, AGE, DIAGNOSIS, assign)
rownames(pheno_) <- NULL
table(pheno_$assign)
# diagnosis data
data <- read_tsv('/datacommons/ochoalab/ssns_gwas/array/patient-data.txt.gz') %>% filter(bristol == TRUE) %>% 
  filter(diagnosis == "SSNS")

outlier_labels = read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/admixture/outlier_labels1.txt")$V1
pheno_data = pheno_ %>% filter(! IID %in% outlier1) %>% select(IID, assign)
table(pheno_data$assign)
# # # merge with Rasheed's list
ssns_pheno_data = merge(pheno_data, data, by.x = "IID", by.y = "id", all.y = TRUE) %>% drop_na(assign)
table(ssns_pheno_data$assign)

# # write each ancestry into txt files
write.table(cbind(0, ssns_pheno_data %>% filter(assign == 1) %>% pull(IID)),
            "/datacommons/ochoalab/ssns_gwas/replication/saige_results/bristol/AF/ssns_ctrl/bristol_ssns_eur.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(cbind(0, ssns_pheno_data %>% filter(assign == 2) %>% pull(IID)),
            "/datacommons/ochoalab/ssns_gwas/replication/saige_results/bristol/AF/ssns_ctrl/bristol_ssns_sas.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(cbind(0, ssns_pheno_data %>% filter(assign == 3) %>% pull(IID)),
            "/datacommons/ochoalab/ssns_gwas/replication/saige_results/bristol/AF/ssns_ctrl/bristol_ssns_afr.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(cbind(0, ssns_pheno_data %>% pull(IID)),
            "/datacommons/ochoalab/ssns_gwas/replication/saige_results/bristol/AF/ssns_ctrl/bristol_ssns_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

# some ref/alt mismatch due to multiallelic cases
# can now run LRT on Bristol cases vs UKBB controls
# also perform LD clump on Bristol genotype data + LRT p-val outputs
## and calculate AC
############
# get bristol af for all snps including imputation bias snps
############
AC_df_clean = control_vs_control %>% 
  filter(!(CHROM %in% list_to_remove$CHROM & POS %in% list_to_remove$POS & A1 %in% list_to_remove$A1 & A2 %in% list_to_remove$A2))

control_snp = AC_df_clean %>% mutate(SNP = paste0(CHROM, ":", POS, ":",  A1, ":", A2)) %>% select(SNP)

write.table(control_snp, paste0(dir, "bristol/AF/", disease_subtype, "/clean_control_snp_id.txt"), 
            row.names=FALSE, quote=FALSE, col.names = FALSE)

# run allele_freq_bristol.q, process outfiles
sas_ac <- read.table(paste0(dir, "bristol/AF/", disease_subtype, "/bristolAC_sas.acount")) 
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
sas_ac_clean = sas_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristol_AC_sas = ALT_FREQS, bristol_AN_sas = OBS_CT)

eur_ac <- read.table(paste0(dir, "bristol/AF/", disease_subtype, "/bristolAC_eur.acount")) 
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristol_AC_eur = ALT_FREQS, bristol_AN_eur = OBS_CT)

afr_ac <- read.table(paste0(dir, "bristol/AF/", disease_subtype, "/bristolAC_afr.acount")) 
colnames(afr_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
afr_ac_clean = afr_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristol_AC_afr = ALT_FREQS, bristol_AN_afr = OBS_CT)

merge1 = merge(sas_ac_clean, eur_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)
merge2 = merge(merge1, afr_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)

write.table(merge2, 
            paste0(dir,  "bristol/AF/", disease_subtype, "/bristol_AC_ssns_allage.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE)

bristol_AC <- read.table(paste0(dir, "bristol/AF/",  disease_subtype, "/bristol_AC_ssns_allage.txt"), header = TRUE) %>% 
  mutate(CHROM = str_remove(CHROM, "chr"))

cases_vs_control = merge(ukbb_AC, bristol_AC, by.x = c("Chromosome", "Position", "A1", "A2"), 
                         by.y = c("CHROM", "POS", "REF", "ALT")) %>% 
  select(Chromosome, Position, A1, A2,  ukbb_AC_afr, ukbb_AN_afr, 
         ukbb_AC_eur, ukbb_AN_eur, ukbb_AC_sas, ukbb_AN_sas, bristol_AC_afr, bristol_AN_afr, 
         bristol_AC_eur, bristol_AN_eur, bristol_AC_sas, bristol_AN_sas) %>% distinct()

# LRT on Bristol vs UKBB (case vs control)
# define function input
# x1/n1: bistrol; x2/n2: UKBB
x1 <- cbind(cases_vs_control$bristol_AC_sas, cases_vs_control$bristol_AC_eur, cases_vs_control$bristol_AC_afr)
n1 <- cbind(cases_vs_control$bristol_AN_sas, cases_vs_control$bristol_AN_eur, cases_vs_control$bristol_AN_afr)
x2 <- cbind(cases_vs_control$ukbb_AC_sas, cases_vs_control$ukbb_AC_eur, cases_vs_control$ukbb_AC_afr)
n2 <- cbind(cases_vs_control$ukbb_AN_sas, cases_vs_control$ukbb_AN_eur, cases_vs_control$ukbb_AN_afr)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals_b <- af_test( x1, n1, x2, n2 )$pval

hist(pvals_b, breaks = 50, main = "Bristol vs UKBB")
qq(pvals_b, main = "Bristol vs UKBB")

output_plot = cbind(cases_vs_control, pvals_b)
write.table(output_plot, paste0(dir, "bristol/AF/annotated_results/replication_pvals_", disease_subtype, ".txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

# define p-value threshold
samplesize = length(pvals_b)
pcut = 0.05/samplesize
cases_vs_control_test = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% 
  dplyr::rename(BP = Position, P = pvals_b, CHR = Chromosome) %>% mutate(SNP = paste0(CHR, ":", BP)) %>% 
  arrange(CHR, BP) %>% mutate(sig = ifelse(P < (pcut), TRUE, FALSE))
table(cases_vs_control_test$sig)

### write file for LD clump (All bristol snps)
bristol_snps = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% dplyr::rename(CHR = Chromosome, POS = Position) %>% 
  mutate(SNP = paste0("chr", CHR, ":", POS, ":", A1, ":", A2)) %>% select(SNP, P = pvals_b)
write.table(bristol_snps, paste0(dir, "bristol/AF/", disease_subtype, "/replication_unclump_bristol.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")
### clump Bristol p-values to get effective number of tests
clump_bristol <- read.table(paste0(dir, "bristol/AF/", disease_subtype, "/clump_bristol_4848.clumped"), header = TRUE) 
nrow(clump_bristol)

# apply new p-value threshold
pval_filter = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < 0.05/nrow(clump_bristol)) %>% 
  dplyr::rename(CHR = Chromosome, POS = Position) %>% 
  mutate(SNP = paste0("chr", CHR, ":", POS, ":", A1, ":", A2))
nrow(pval_filter)
length(intersect(clump_bristol$SNP, pval_filter$SNP))

# annotations from discovery
anno_rsid = read.csv(paste0("/datacommons/ochoalab/ssns_gwas/saige/annotated/", disease_subtype, ".csv")) %>% 
  dplyr::rename(original_PVAL = PVAL) %>% select(-BETA, -SE, -Tstat, -VAR, -AF_case, -N_case, -AF_ctrl, -N_ctrl, -X)
sumstat_filter = anno_rsid %>% filter(Chromosome %in% pval_filter$CHR & Position %in% pval_filter$POS) 
replication_merge = merge(sumstat_filter, pval_filter %>% 
                            select(CHR, POS, A1, A2, Bristol_LRT_PVAL = pvals_b), 
                          by.x = c("Chromosome", "Position", "A1", "A2"), by.y = c("CHR", "POS", "A1", "A2"))

write.table(replication_merge, paste0(dir, "bristol/AF/annotated_results/bristol_replication_", disease_subtype, ".txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

# annotate snps according to ld clumped number of independent snps

cases_vs_control_pval = cases_vs_control_test %>% mutate(sig = ifelse(P < (0.05/nrow(clump_bristol)), TRUE, FALSE))
table(cases_vs_control_pval$sig)
# rsid_list = c("rs9273807", "rs2856695", "rs114032596", "rs567166387", "rs11860603", 
#               "rs184934338", "rs150918987", "rs2855812")

snp_list = c("6:32652506", "10:28810849", "16:11077745", "6:32689478", "6:31361670")

snp_list = c("17:56125448", "8:109898199", "8:107771874", "8:135891550", "7:3014907", "15:101859143",
             "8:142322618", "10:45904808", "1:218141640")

ssns_repli_table = cases_vs_control_test %>% filter(SNP %in% snp_list)

