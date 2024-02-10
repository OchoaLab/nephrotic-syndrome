library(tidyverse)
library(genio)
library(qqman)
dir = '/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/'
setwd( dir )

#############################
disease_subtype = "ns_ctrl"
ancestry = "eur"
#dir = '/datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/'
#############################

# step 1: extract gnomad SNP id's
gnomad_AC <- read.table(paste0(dir, "gnomad/outfile_new/gnomad-genome_", disease_subtype, "_", ancestry, ".txt"), sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE) %>% 
  filter(AC_nfe != 0)
gnomad_snp = gnomad_AC %>% mutate(SNP = paste0(CHROM, ":", POS, ":",  REF, ":", ALT)) %>% select(SNP)
write.table(gnomad_snp, paste0(dir, "gnomad/outfile_new/snp-id/", disease_subtype, "_", ancestry, "_snp_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)

# step 2: run allele_freq_discovery.q on NS data to remove false positive SNPs
# load phenotype and fixed covariates file
# write id files to calculate AC per ancestry
# data <- read_tsv( '/datacommons/ochoalab/ssns_gwas/imputed/patient-data.txt.gz', col_types = 'cccccdciiii' )
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

# step 3: run AC_processing.q (discovery_AC_processing.R)
eur_ac <- read.table(paste0(dir, disease_subtype, "/", ancestry, "/discoveryAC_eur_gnomadsnps.acount")) 
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(disc_AC_eur = ALT_FREQS, disc_AN_eur = OBS_CT)

write.table(eur_ac_clean, paste0(dir, disease_subtype, "/", ancestry, "/discovery_AC_eur.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE)

# step 4: Merge discovery AC and gnomad AC
discovery_AC <- read.table(paste0(dir, disease_subtype, "/", ancestry, "/discovery_AC_eur.txt"), 
                           header = TRUE, sep = " ")
AC_df = merge(gnomad_AC, discovery_AC, by = c("CHROM", "POS", "REF", "ALT"))

# step 5: AF test for discovery vs gnomad
library(popgeninfer)
# define function input
# x1/n1: discovery set; x2/n2: gnomad
x1 <- AC_df$disc_AC_eur
n1 <- AC_df$disc_AN_eur
x2 <- AC_df$AC_nfe
n2 <- AC_df$AN_nfe
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals <- af_test( x1, n1, x2, n2 )$pval
# define p-value threshold
samplesize = length(pvals)
pcut = 0.05/samplesize
list_to_remove = cbind(AC_df, pvals) %>% arrange(pvals) %>% filter(pvals < pcut) 
list_to_remove

hist(pvals, breaks = 50, main = "p-vals likelihood ratio test (all ancestry)")
qq(pvals, main = "QQ-plot for all ancestry")

# step 6: Bristol AF test with gnomad with clean set of SNPs: 
## write SNPs into file 
## write per ancestry id's
# library(rio)
# # reading data from all sheets
# data <- import_list("../PHENOTYPE_BRISTOL_REPLICATION PLANS.xlsx") 
# NS_sheet = data$ALL[,1:5]
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
# ns_pheno_data = merge(pheno_data, NS_sheet, by.x = "IID", by.y = "AcquisitionNumber", all.y = TRUE) %>% drop_na(assign)
# table(ns_pheno_data$assign)
# # write each ancestry into txt files
# write.table(cbind(0, ns_pheno_data %>% filter(assign == 1) %>% pull(IID)),
#             "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ns_ctr_all/bristol_ns_euro_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, ns_pheno_data %>% filter(assign == 2) %>% pull(IID)),
#             "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ns_ctr_all/bristol_ns_sas_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, ns_pheno_data %>% filter(assign == 3) %>% pull(IID)),
#             "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ns_ctr_all/bristol_ns_af_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)
# write.table(cbind(0, ns_pheno_data %>% pull(IID)),
#             "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ns_ctr_all/bristol_ns_allage.txt", row.names=FALSE, quote=FALSE, col.names = FALSE)

# some ref/alt mismatch due to multiallelic cases
# can now run LRT on Bristol cases vs Gnomad controls
# also perform LD clump on Bristol genotype data + LRT p-val outputs
## and calculate AC
AC_df_clean = AC_df %>% 
  filter(!(CHROM %in% list_to_remove$CHROM & POS %in% list_to_remove$POS & REF %in% list_to_remove$REF & ALT %in% list_to_remove$ALT))

control_snp = AC_df_clean %>% mutate(SNP = paste0(CHROM, ":", POS, ":",  REF, ":", ALT)) %>% select(SNP)
write.table(control_snp,
            paste0(dir, disease_subtype,  "/", ancestry, "/control_snp_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)
# run allele_freq_bristol.q, process outfiles
eur_ac <- read.table(paste0(dir, disease_subtype,  "/", ancestry, "/bristolAC_eur.acount")) 
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHROM", "POS", "A1", "A2")) %>% select(-A1, -A2) %>% dplyr::rename(bristol_AC_eur = ALT_FREQS, bristol_AN_eur = OBS_CT)

write.table(eur_ac_clean, 
            paste0(dir, disease_subtype,  "/", ancestry, "/bristol_AC_ns_eur.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE)

bristol_AC <- read.table(paste0(dir, disease_subtype,  "/", ancestry,  "/bristol_AC_ns_eur.txt"), header = TRUE) 
cases_vs_control = merge(gnomad_AC, bristol_AC, by = c("CHROM", "POS", "REF", "ALT"))

# LRT on Bristol vs Gnomad (case vs control)
# define function input
# x1/n1: bistrol; x2/n2: gnomad
x1 <- cases_vs_control$bristol_AC_eur
n1 <- cases_vs_control$bristol_AN_eur
x2 <- cases_vs_control$AC_nfe
n2 <- cases_vs_control$AN_nfe
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals_b <- af_test( x1, n1, x2, n2 )$pval
# define p-value threshold
samplesize = length(pvals_b)
pcut = 0.05/samplesize

# define p-value threshold
output_plot = cbind(cases_vs_control, pvals_b)
write.table(output_plot, paste0(dir, "/annotated_results/replication_pvals_", disease_subtype, "_", ancestry, ".txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

unclumped_sig = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < pcut) 
nrow(unclumped_sig)

hist(pvals_b, breaks = 50, main = "p-vals likelihood ratio test (all ancestry)")
qq(pvals_b, main = "QQ-plot for all ancestry")

### write file for LD clump (All bristol snps)
bristol_snps = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>%
  mutate(SNP = paste0(CHROM, ":", POS, ":", REF, ":", ALT)) %>% select(SNP, P = pvals_b)
write.table(bristol_snps, paste0(dir, disease_subtype, "/", ancestry, "/replication_unclump_bristol.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")
### clump Bristol p-values to get effective number of tests
clump_bristol <- read.table(paste0(dir, disease_subtype, "/", ancestry, "/clump_bristol_1440.clumped"), header = TRUE) 
nrow(clump_bristol)

# apply new p-value threshold
pval_filter = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < 0.05/nrow(clump_bristol)) %>% 
  mutate(SNP = paste0(CHROM, ":", POS, ":", REF, ":", ALT))
nrow(pval_filter)
length(intersect(clump_bristol$SNP, pval_filter$SNP))
intersect_list = intersect(clump_bristol$SNP, pval_filter$SNP)
independent_regions = pval_filter %>% filter(SNP %in% intersect_list) %>% select(CHROM, POS, REF, ALT, ID, SNP, pvals_b) %>% 
  mutate(CHROM = str_remove(CHROM, "chr"))

# write results to SNPNexus
SNP_nexus = independent_regions %>% mutate(type = "Chromosome", Strand = 1) %>% 
  select(type, Id = CHROM, Position = POS, Allele1 = REF, Allele2 = ALT, Strand) %>% distinct()

write.table(SNP_nexus, paste0(dir, disease_subtype, "/", ancestry, "/bristol_snp_nexus.txt"), 
            row.names=FALSE, quote=FALSE, col.names = FALSE, sep = "\t")

neargens = read.table(paste0(dir, disease_subtype, "/", ancestry, "/near_gens_bristol_", disease_subtype, "_", ancestry, ".txt"), header = TRUE, sep = "\t") %>% select(-Chromosome, -Position)
gencoord = read.table(paste0(dir, disease_subtype, "/", ancestry, "/gen_coords_bristol_", disease_subtype, "_", ancestry, ".txt"), header = TRUE, sep = "\t") %>% 
  select(Variation.ID, Chromosome, Position, A1 = REF.Allele, A2 = ALT.Allele..IUPAC., rsid = dbSNP)
merge_nexus = merge(gencoord, neargens, by = "Variation.ID") %>% select(-Variation.ID)

# include summary stats in output (original gmmat)
#ancestry_all = read.table(paste0("/datacommons/ochoalab/ssns_gwas/imputed/", disease_subtype, "/mac20-glmm-score.txt" ),  sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
ancestry_all = read.table(paste0("/datacommons/ochoalab/ssns_gwas/imputed/",  ancestry, "/mac20-glmm-score.txt" ),  sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)

sumstat_filter = ancestry_all %>% filter(CHR %in% independent_regions$CHROM & POS %in% independent_regions$POS) %>% 
  select(CHR, BP = POS, A1 = REF, A2 = ALT, original_AF = AF, original_SCORE = SCORE, original_VAR = VAR, original_PVAL = PVAL)
ind_reg_merge = merge(independent_regions %>% dplyr::rename(CHR = CHROM, BP = POS, A1 = REF, A2 = ALT), 
                      sumstat_filter, by = c("CHR", "BP", "A1", "A2")) %>% select(Chromosome = CHR, Position = BP, A1, A2, original_AF, original_SCORE, 
         original_VAR, original_PVAL, Bristol_LRT_PVAL = pvals_b)

ind_reg_anno = merge(ind_reg_merge, merge_nexus, by = c("Chromosome", "Position", "A1", "A2")) %>% 
  select(Chromosome, Position, A1, A2, rsid, everything())

write.table(ind_reg_anno, paste0(dir, "/annotated_results/bristol_replication_", disease_subtype, "_", ancestry, ".txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

independent_regions %>% filter(ID == "rs9274564")


# annotate snps according to ld clumped number of independent snps

cases_vs_control_pval = cbind(cases_vs_control, pvals_b) %>% mutate(sig = ifelse(pvals_b < (0.05/nrow(clump_bristol)), TRUE, FALSE))
table(cases_vs_control_pval$sig)


SNP_nexus = cases_vs_control_pval %>% filter(sig == TRUE) %>% select(CHR = CHROM, BP = POS, REF, ALT, ID, P = pvals_b) %>% 
  mutate(type = "Chromosome", Strand = 1, CHR = str_remove(CHR, "chr")) %>% 
  select(type, Id = CHR, Position = BP, Allele1 = REF, Allele2 = ALT, Strand) %>% distinct()
write.table(SNP_nexus, paste0(dir, disease_subtype, "/bristol_snp_nexus_allsig_eur.txt"), 
            row.names=FALSE, quote=FALSE, col.names = FALSE, sep = "\t")

neargens = read.table(paste0(dir, disease_subtype, "/", ancestry, "/near_gens_bristol_", disease_subtype, "_", ancestry, "_allsig.txt"), header = TRUE, sep = "\t") %>% select(-Chromosome, -Position)
gencoord = read.table(paste0(dir, disease_subtype, "/", ancestry,  "/gen_coords_bristol_", disease_subtype, "_", ancestry, "_allsig.txt"), header = TRUE, sep = "\t") %>% 
  select(Variation.ID, Chromosome, Position, A1 = REF.Allele, A2 = ALT.Allele..IUPAC., rsid = dbSNP)
merge_nexus = merge(gencoord, neargens, by = "Variation.ID") %>% select(-Variation.ID)

# data from discovery set
sumstat_filter = ancestry_all %>% filter(CHR %in% SNP_nexus$Id & POS %in% SNP_nexus$Position) %>% 
  select(CHR, BP = POS, A1 = REF, A2 = ALT, original_AF = AF, original_SCORE = SCORE, original_VAR = VAR, original_PVAL = PVAL)
all_snps = cases_vs_control_pval %>% filter(sig == TRUE) %>% select(CHR = CHROM, BP = POS, A1 = REF, A2 = ALT, rsid = ID, 
                                                                    Bristol_LRT_PVAL = pvals_b) %>% mutate(CHR = str_remove(CHR, "chr"))
allsnps_merge = merge(all_snps, sumstat_filter, 
                      by = c("CHR", "BP", "A1", "A2")) %>% select(Chromosome = CHR, Position = BP, A1, A2, original_AF, original_SCORE, 
                                                                  original_VAR, original_PVAL, Bristol_LRT_PVAL)

allsnps_merge_anno = merge(allsnps_merge, merge_nexus, by = c("Chromosome", "Position", "A1", "A2")) %>% 
  select(Chromosome, Position, A1, A2, rsid, everything()) %>% arrange(original_PVAL)

write.table(allsnps_merge_anno, paste0(dir, "/annotated_results/bristol_replication_", disease_subtype, "_", ancestry,  "_sub_allsigsnps.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

