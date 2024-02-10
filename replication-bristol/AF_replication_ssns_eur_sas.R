# Bristol replication is only performed on EUR and SAS ancestry due to low sample size in AFR
library(tidyverse)
library(genio)
library(qqman)
dir = '/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/'
setwd( dir )

#############################
disease_subtype = "ssns_ctrl"
ancestry = "all"
#dir = '/datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/'
#############################

# step 1: extract gnomad SNP id's
gnomad_AC <- read.table(paste0(dir, "gnomad/outfile_new/gnomad-genome_", disease_subtype, "_", ancestry, ".txt"), sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
gnomad_snp = gnomad_AC %>% mutate(SNP = paste0(CHROM, ":", POS, ":",  REF, ":", ALT)) %>% select(SNP)
write.table(gnomad_snp, paste0(dir, "gnomad/outfile_new/snp-id/", disease_subtype, "_", ancestry, "_snp_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)

# step 2: run allele_freq_discovery.q on NS data to remove false positive SNPs
# load phenotype and fixed covariates file
# write id files to calculate AC per ancestry
data <- read_tsv( '/datacommons/ochoalab/ssns_gwas/imputed/patient-data.txt.gz', col_types = 'cccccdciiii' )
sas_id = data %>% filter(diagnosis == "Control" & ancestry == "sas") %>% select(id)
eur_id = data %>% filter(diagnosis == "Control" & ancestry == "eur") %>% select(id)
all_id = data %>% filter(diagnosis == "Control") %>% select(id)

write.table(cbind(0, sas_id), paste0(dir, disease_subtype, "/ctrl_sas_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(cbind(0, eur_id), paste0(dir, disease_subtype, "/ctrl_eur_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)
write.table(cbind(0, all_id), paste0(dir, disease_subtype, "/ctrl_all_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)

# step 3: run AC_processing.q (discovery_AC_processing.R)
sas_ac <- read.table(paste0(dir, disease_subtype, "/discoveryAC_sas_gnomadsnps.acount")) 
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
sas_ac_clean = sas_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(disc_AC_sas = ALT_FREQS, disc_AN_sas = OBS_CT)

eur_ac <- read.table(paste0(dir, disease_subtype, "/discoveryAC_eur_gnomadsnps.acount")) 
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(disc_AC_eur = ALT_FREQS, disc_AN_eur = OBS_CT)

merge1 = merge(sas_ac_clean, eur_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)

write.table(merge1, paste0(dir, disease_subtype, "/discovery_AC_eur_sas.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE)

# step 4: Merge discovery AC and gnomad AC
discovery_AC <- read.table(paste0(dir, disease_subtype, "/discovery_AC_eur_sas.txt"), 
                           header = TRUE, sep = " ")
AC_df = merge(gnomad_AC, discovery_AC, by = c("CHROM", "POS", "REF", "ALT"))

# step 5: AF test for discovery vs gnomad
library(popgeninfer)
# define function input
# x1/n1: discovery set; x2/n2: gnomad
x1 <- cbind(AC_df$disc_AC_sas, AC_df$disc_AC_eur)
n1 <- cbind(AC_df$disc_AN_sas, AC_df$disc_AN_eur)
x2 <- cbind(AC_df$AC_sas, AC_df$AC_nfe)
n2 <- cbind(AC_df$AN_sas, AC_df$AN_nfe)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals <- af_test( x1, n1, x2, n2 )$pval
# sort by p-values
control_test = cbind(AC_df, pvals) %>% arrange(pvals) %>% 
  dplyr::rename(CHR = CHROM, BP = POS, P = pvals) %>% mutate(CHR = as.numeric(str_remove(CHR, "chr")), SNP = paste0(CHR, ":", BP)) %>% 
  arrange(CHR, BP) %>% mutate(sig = ifelse(P < pcut, TRUE, FALSE))
table(control_test$sig)
# define p-value threshold
samplesize = length(pvals)
pcut = 0.05/samplesize
list_to_remove = cbind(AC_df, pvals) %>% arrange(pvals) %>% filter(pvals < pcut) 
list_to_remove

hist(pvals, breaks = 50, main = "discovery control vs gnomad control (eur and sas)")
qq(pvals, main = "discovery control vs gnomad control")

# step 6: Bristol AF test with gnomad with clean set of SNPs: 
## write SNPs into file 
## write per ancestry id's
# library(rio)
# # reading data from all sheets
# data <- import_list("../PHENOTYPE_BRISTOL_REPLICATION PLANS.xlsx") 
# SSNS_sheet = data$SSNS[,1:5]
# SSNS_age = data$ALL[,1:5] %>% mutate(AGE = as.numeric(AGE)) %>%  filter( DIAGNOSIS == "SSNS")
# 
# # write subset id's with phenotype
# pheno_ <- read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/admixture/pca1_pheno_df.txt", sep = "\t", header = TRUE) %>% 
#   select(IID, Sex, RACE, AGE, DIAGNOSIS, assign)
# rownames(pheno_) <- NULL
# table(pheno_$assign)
# outlier1 = read.table("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/admixture/outlier_labels1.txt")$V1
# pheno_data = pheno_ %>% filter(! IID %in% outlier1) %>% select(IID, assign)
# table(pheno_data$assign)
# # merge with Rasheed's list
# ssns_pheno_data = merge(pheno_data, SSNS_age, by.x = "IID", by.y = "AcquisitionNumber", all.y = TRUE) %>% drop_na(assign)
# table(ssns_pheno_data$assign)
# # write each ancestry into txt files
# write.table(cbind(0, ssns_pheno_data %>% filter(assign == 1) %>% pull(IID)),
#             "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/bristol_ssns_euro_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, ssns_pheno_data %>% filter(assign == 2) %>% pull(IID)),
#             "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/bristol_ssns_sas_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, ssns_pheno_data %>% filter(assign == 3) %>% pull(IID)),
#             "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/bristol_ssns_af_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, ssns_pheno_data %>% pull(IID)),
#             "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/bristol_ssns_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

# some ref/alt mismatch due to multiallelic cases
# can now run LRT on Bristol cases vs Gnomad controls
# also perform LD clump on Bristol genotype data + LRT p-val outputs
## and calculate AC
AC_df_clean = AC_df %>% 
  filter(!(CHROM %in% list_to_remove$CHROM & POS %in% list_to_remove$POS & REF %in% list_to_remove$REF & ALT %in% list_to_remove$ALT))

control_snp = AC_df_clean %>% mutate(SNP = paste0(CHROM, ":", POS, ":",  REF, ":", ALT)) %>% select(SNP)
write.table(control_snp,
            "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctrl/control_snp_id_eur_sas.txt", 
            row.names=FALSE, quote=FALSE, col.names = FALSE)
# run allele_freq_bristol.q, process outfiles
sas_ac <- read.table(paste0(dir, disease_subtype, "/bristolAC_sas_all_new.acount")) 
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
sas_ac_clean = sas_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristol_AC_sas = ALT_FREQS, bristol_AN_sas = OBS_CT)

eur_ac <- read.table(paste0(dir, disease_subtype, "/bristolAC_eur_all_new.acount")) 
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristol_AC_eur = ALT_FREQS, bristol_AN_eur = OBS_CT)

merge1 = merge(sas_ac_clean, eur_ac_clean, by = c("CHROM", "POS", "REF", "ALT"), all = TRUE)

write.table(merge1, 
            paste0(dir, disease_subtype, "/bristol_AC_ssns_allage_eur_sas.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE)

bristol_AC <- read.table(paste0(dir, disease_subtype, "/bristol_AC_ssns_allage_eur_sas.txt"), header = TRUE) 
cases_vs_control = merge(gnomad_AC, bristol_AC, by = c("CHROM", "POS", "REF", "ALT"))

# LRT on Bristol vs Gnomad (case vs control)
# define function input
# x1/n1: bistrol; x2/n2: gnomad
x1 <- cbind(cases_vs_control$bristol_AC_sas, cases_vs_control$bristol_AC_eur)
n1 <- cbind(cases_vs_control$bristol_AN_sas, cases_vs_control$bristol_AN_eur)
x2 <- cbind(cases_vs_control$AC_sas, cases_vs_control$AC_nfe)
n2 <- cbind(cases_vs_control$AN_sas, cases_vs_control$AN_nfe)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals_b <- af_test( x1, n1, x2, n2 )$pval
# define p-value threshold
samplesize = length(pvals_b)
pcut = 0.05/samplesize
cases_vs_control_test = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% 
  dplyr::rename(CHR = CHROM, BP = POS, P = pvals_b) %>% mutate(CHR = as.numeric(str_remove(CHR, "chr")), SNP = paste0(CHR, ":", BP)) %>% 
  arrange(CHR, BP) %>% mutate(sig = ifelse(P < pcut, TRUE, FALSE))
table(cases_vs_control_test$sig)


# define p-value threshold
# list_to_remove = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < pcut) 
# list_to_remove

hist(pvals_b, breaks = 50, main = "Bristol vs Gnomad (eur and sas)")
qq(pvals_b, main = "Bristol vs Gnomad")

output_plot = cbind(cases_vs_control, pvals_b)
write.table(output_plot, paste0(dir, "/annotated_results/replication_pvals_", disease_subtype, "_", ancestry, "_eur_sas.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

### write file for LD clump (All bristol snps)
bristol_snps = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>%
  mutate(SNP = paste0(CHROM, ":", POS, ":", REF, ":", ALT)) %>% select(SNP, P = pvals_b)
write.table(bristol_snps, paste0(dir, disease_subtype, "/replication_unclump_bristol_allage_eur_sas.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")
### clump Bristol p-values to get effective number of tests
clump_bristol <- read.table(paste0(dir, disease_subtype, "/clump_bristol_6350.clumped"), header = TRUE) 
nrow(clump_bristol)

# apply new p-value threshold
pval_filter = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < 0.05/nrow(clump_bristol)) %>% 
  mutate(SNP = paste0(CHROM, ":", POS, ":", REF, ":", ALT))
length(intersect(clump_bristol$SNP, pval_filter$SNP))
intersect_list = intersect(clump_bristol$SNP, pval_filter$SNP)
independent_regions = pval_filter %>% filter(SNP %in% intersect_list) %>% select(CHROM, POS, REF, ALT, ID, SNP, pvals_b) %>% 
  mutate(CHROM = str_remove(CHROM, "chr"))

# write results to SNPNexus
SNP_nexus = independent_regions %>% mutate(type = "Chromosome", Strand = 1) %>% 
  select(type, Id = CHROM, Position = POS, Allele1 = REF, Allele2 = ALT, Strand) %>% distinct()

write.table(SNP_nexus, paste0(dir, disease_subtype, "/bristol_snp_nexus_new.txt"), 
            row.names=FALSE, quote=FALSE, col.names = FALSE, sep = "\t")

neargens = read.table(paste0(dir, disease_subtype, "/near_gens_bristol_", disease_subtype, "_", ancestry, "_eur_sas.txt"), header = TRUE, sep = "\t") %>% select(-Chromosome, -Position)
gencoord = read.table(paste0(dir, disease_subtype, "/gen_coords_bristol_", disease_subtype, "_", ancestry, "_eur_sas.txt"), header = TRUE, sep = "\t") %>% 
  select(Variation.ID, Chromosome, Position, A1 = REF.Allele, A2 = ALT.Allele..IUPAC., rsid = dbSNP)
merge_nexus = merge(gencoord, neargens, by = "Variation.ID") %>% select(-Variation.ID)

# include summary stats in output (original gmmat)
ancestry_all = read.table(paste0("/datacommons/ochoalab/ssns_gwas/imputed/", disease_subtype, "/mac20-glmm-score.txt" ),  sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
sumstat_filter = ancestry_all %>% filter(CHR %in% independent_regions$CHROM & POS %in% independent_regions$POS) %>% 
  select(CHR, BP = POS, A1 = REF, A2 = ALT, original_AF = AF, original_SCORE = SCORE, original_VAR = VAR, original_PVAL = PVAL)
ind_reg_merge = merge(independent_regions %>% dplyr::rename(CHR = CHROM, BP = POS, A1 = REF, A2 = ALT), 
                      sumstat_filter, by = c("CHR", "BP", "A1", "A2")) %>% select(Chromosome = CHR, Position = BP, A1, A2, original_AF, original_SCORE, 
         original_VAR, original_PVAL, Bristol_LRT_PVAL = pvals_b) %>% arrange(original_PVAL)

ind_reg_anno = merge(ind_reg_merge, merge_nexus, by = c("Chromosome", "Position", "A1", "A2")) %>% 
  select(Chromosome, Position, A1, A2, rsid, everything())

write.table(ind_reg_anno, paste0(dir, "/annotated_results/bristol_replication_", disease_subtype, "_eur_sas.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

write.table(cases_vs_control_test %>% select(CHR, BP, REF, ALT, ID, P), paste0(dir, "/annotated_results/replication_pvals_", disease_subtype, "_", ancestry, ".txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

# annotate snps according to ld clumped number of independent snps

cases_vs_control_pval = cases_vs_control_test %>% mutate(sig = ifelse(P < (0.05/nrow(clump_bristol)), TRUE, FALSE))
table(cases_vs_control_pval$sig)
rsid_list = c("rs3844313", "rs2856695", "rs114032596", "rs12925642", "rs11213573", 
              "rs62072970", "rs139087832", "rs151109286", "rs140718574", "rs150401529")

ssns_repli_table = cases_vs_control_test %>% filter(ID %in% rsid_list)

SNP_nexus = cases_vs_control_pval %>% filter(sig == TRUE) %>% select(CHR, BP, REF, ALT, ID, SNP, P) %>% 
  mutate(type = "Chromosome", Strand = 1) %>% 
  select(type, Id = CHR, Position = BP, Allele1 = REF, Allele2 = ALT, Strand) %>% distinct()
write.table(SNP_nexus, paste0(dir, disease_subtype, "/bristol_snp_nexus_allsig.txt"), 
            row.names=FALSE, quote=FALSE, col.names = FALSE, sep = "\t")

neargens = read.table(paste0(dir, disease_subtype, "/near_gens_bristol_", disease_subtype, "_", ancestry, "_eur_sas_allsig.txt"), header = TRUE, sep = "\t") %>% select(-Chromosome, -Position)
gencoord = read.table(paste0(dir, disease_subtype, "/gen_coords_bristol_", disease_subtype, "_", ancestry, "_eur_sas_allsig.txt"), header = TRUE, sep = "\t") %>% 
  select(Variation.ID, Chromosome, Position, A1 = REF.Allele, A2 = ALT.Allele..IUPAC., rsid = dbSNP)
merge_nexus = merge(gencoord, neargens, by = "Variation.ID") %>% select(-Variation.ID)

# data from discovery set
sumstat_filter = ancestry_all %>% filter(CHR %in% SNP_nexus$Id & POS %in% SNP_nexus$Position) %>% 
  select(CHR, BP = POS, A1 = REF, A2 = ALT, original_AF = AF, original_SCORE = SCORE, original_VAR = VAR, original_PVAL = PVAL)
all_snps = cases_vs_control_pval %>% filter(sig == TRUE) %>% select(CHR, BP, A1 = REF, A2 = ALT, rsid = ID, 
                                                                    Bristol_LRT_PVAL = P)
allsnps_merge = merge(all_snps, sumstat_filter, 
                      by = c("CHR", "BP", "A1", "A2")) %>% select(Chromosome = CHR, Position = BP, A1, A2, original_AF, original_SCORE, 
                                                                  original_VAR, original_PVAL, Bristol_LRT_PVAL)

allsnps_merge_anno = merge(allsnps_merge, merge_nexus, by = c("Chromosome", "Position", "A1", "A2")) %>% 
  select(Chromosome, Position, A1, A2, rsid, everything()) %>% arrange(original_PVAL)

write.table(allsnps_merge_anno, paste0(dir, "/annotated_results/bristol_replication_", disease_subtype, "_eur_sas_allsigsnps.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

