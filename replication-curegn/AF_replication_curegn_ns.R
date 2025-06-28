library(tidyverse)
library(genio)
library(qqman)
library(popgeninfer)
dir = '/datacommons/ochoalab/ssns_gwas/replication/saige_results/'
setwd( dir )

#############################
disease_subtype = "ns_ctrl"
ancestry = "all"

# step 1: extract gnomad SNP id's
gnomad_AC <- read.table(paste0(dir, "gnomad/outfile/gnomad-genome_", disease_subtype, "_clean.txt"), sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
gnomad_snp = gnomad_AC %>% mutate(SNP = paste0(CHROM, ":", POS)) %>% mutate(SNP = str_remove(SNP, "chr")) %>% select(SNP)
write.table(gnomad_snp, paste0(dir, "curegn/AF/ns/gnomad_snp_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)

# curegn only has cases, no controls, so skip step 2 of removing imputation biased snps from control vs control
# run LRT on CureGN cases vs Gnomad controls
# run allele_freq_curegn.q, process outfiles
sas_ac <- read.table(paste0(dir, "curegn/AF/ns/curegnAC_sas.acount")) 
colnames(sas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT") 
sas_ac_clean = sas_ac %>% 
  separate(ID, c("CHR", "POS")) %>% select(-CHROM) %>% dplyr::rename(curegn_AC_sas = ALT_FREQS, curegn_AN_sas = OBS_CT)

eur_ac <- read.table(paste0(dir, "curegn/AF/ns/curegnAC_eur.acount")) 
colnames(eur_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eur_ac_clean = eur_ac  %>%
  separate(ID, c("CHR", "POS")) %>% select(-CHROM)  %>% dplyr::rename(curegn_AC_eur = ALT_FREQS, curegn_AN_eur = OBS_CT)

afr_ac <- read.table(paste0(dir, "curegn/AF/ns/curegnAC_afr.acount")) 
colnames(afr_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
afr_ac_clean = afr_ac  %>%
  separate(ID, c("CHR", "POS")) %>% select(-CHROM)  %>% dplyr::rename(curegn_AC_afr = ALT_FREQS, curegn_AN_afr = OBS_CT)

eas_ac <- read.table(paste0(dir, "curegn/AF/ns/curegnAC_eas.acount")) 
colnames(eas_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
eas_ac_clean = eas_ac  %>%
  separate(ID, c("CHR", "POS")) %>% select(-CHROM)  %>% dplyr::rename(curegn_AC_eas = ALT_FREQS, curegn_AN_eas = OBS_CT)


merge1 = merge(sas_ac_clean, eur_ac_clean, by = c("CHR", "POS", "REF", "ALT"), all = TRUE)
merge2 = merge(merge1, afr_ac_clean, by = c("CHR", "POS", "REF", "ALT"), all = TRUE)
merge3 = merge(merge2, eas_ac_clean, by = c("CHR", "POS", "REF", "ALT"), all = TRUE)

write.table(merge3, paste0(dir, "curegn/AF/ns/curegn_AC_ns.txt"),
            row.names=FALSE, quote=FALSE, col.names = TRUE)

curegn_AC <- read.table(paste0(dir, "curegn/AF/ns/curegn_AC_ns.txt"), header = TRUE) 
cases_vs_control = merge(gnomad_AC, curegn_AC %>% mutate(CHROM = paste0("chr", CHR)), by = c("CHROM", "POS", "REF", "ALT")) 

# LRT on curegn vs Gnomad (case vs control)
# define function input
# x1/n1: curegn; x2/n2: gnomad
x1 <- cbind(cases_vs_control$curegn_AC_sas, cases_vs_control$curegn_AC_eur, cases_vs_control$curegn_AC_afr, cases_vs_control$curegn_AC_eas)
n1 <- cbind(cases_vs_control$curegn_AN_sas, cases_vs_control$curegn_AN_eur, cases_vs_control$curegn_AN_afr, cases_vs_control$curegn_AN_eas )
x2 <- cbind(cases_vs_control$AC_sas, cases_vs_control$AC_nfe, cases_vs_control$AC_afr, cases_vs_control$AC_eas)
n2 <- cbind(cases_vs_control$AN_sas, cases_vs_control$AN_nfe, cases_vs_control$AN_afr, cases_vs_control$AN_eas)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals_b <- af_test( x1, n1, x2, n2 )$pval
# define p-value threshold
samplesize = length(pvals_b)
pcut = 0.05/samplesize
cases_vs_control_test = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% select(-CHR) %>% 
  dplyr::rename(CHR = CHROM, BP = POS, P = pvals_b) %>% #mutate(CHR = as.numeric(str_remove(CHR, "chr")), SNP = paste0(CHR, ":", BP)) %>% 
  arrange(CHR, BP) %>% mutate(sig = ifelse(P < pcut, TRUE, FALSE))

output_plot = cbind(cases_vs_control, pvals_b)
write.table(output_plot, paste0(dir, "curegn/AF/annotated_results/replication_pvals_", disease_subtype,".txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

# define p-value threshold
unclumped_sig = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < pcut) 
nrow(unclumped_sig)

hist(pvals_b, breaks = 50, main = "curegn NS vs gnomad controls")
qq(pvals_b, main = "curegn NS vs gnomad controls")

### write file for LD clump
curegn_snps = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% mutate(CHROM = str_remove(CHROM, "chr")) %>% 
  mutate(SNP = paste0(CHROM, ":", POS)) %>% select(SNP, P = pvals_b)
write.table(curegn_snps, paste0(dir, "curegn/AF/ns/replication_unclump_curegn_snps.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")
### clump p-values to get effective number of tests
clump_curegn <- read.table(paste0(dir,  "curegn/AF/ns/clump_curegn_5746.clumped"), header = TRUE) 
nrow(clump_curegn)

# apply new p-value threshold
pval_filter = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < 0.05/nrow(clump_curegn)) %>% 
  mutate(SNP = paste0(CHROM, ":", POS, ":", REF, ":", ALT))
nrow(pval_filter)
length(intersect(clump_curegn$SNP, pval_filter$SNP))
intersect_list = intersect(clump_curegn$SNP, pval_filter$SNP)
independent_regions = pval_filter %>% filter(SNP %in% intersect_list) %>% select(CHROM, POS, REF, ALT, ID, SNP, pvals_b) %>% 
  mutate(CHROM = str_remove(CHROM, "chr"))


# include summary stats in output (original gmmat)
ancestry_all = read.table("/datacommons/ochoalab/ssns_gwas/saige/ns_ctrl/saige_output.txt" ,  sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# annotate snps according to ld clumped number of independent snps

cases_vs_control_pval = cases_vs_control_test %>% mutate(sig = ifelse(P < (0.05/nrow(clump_curegn)), TRUE, FALSE))
table(cases_vs_control_pval$sig)
rsid_list = c("rs17843604", "rs114032596", "rs11860603", "rs567166387", "rs3129713", "rs2855812")

ns_repli_table = cases_vs_control_test %>% filter(ID %in% rsid_list) %>% select(ID, P, sig) 

SNP_nexus = cases_vs_control_pval %>% filter(sig == TRUE) %>% select(CHR, BP, REF, ALT, ID, P) %>% 
  mutate(CHR = str_remove(CHR, "chr")) %>% 
  mutate(type = "Chromosome", Strand = 1) %>% 
  select(type, Id = CHR, Position = BP, Allele1 = REF, Allele2 = ALT, Strand) %>% distinct()
write.table(SNP_nexus, paste0(dir,"curegn/AF/ns/curegn_snp_nexus_allsig.txt"), 
            row.names=FALSE, quote=FALSE, col.names = FALSE, sep = "\t")

neargens = read.table(paste0(dir,"curegn/AF/ns/near_gens.txt"), header = TRUE, sep = "\t") %>% select(-Chromosome, -Position)
gencoord = read.table(paste0(dir, "curegn/AF/ns/gen_coords.txt"), header = TRUE, sep = "\t") %>% 
  select(Variation.ID, Chromosome, Position, A1 = REF.Allele, A2 = ALT.Allele..IUPAC., rsid = dbSNP)
merge_nexus = merge(gencoord, neargens, by = "Variation.ID") %>% select(-Variation.ID)

# data from discovery set
sumstat_filter = ancestry_all %>% filter(CHR %in% SNP_nexus$Id & POS %in% SNP_nexus$Position) %>% 
  select(Chromosome = CHR, Position = POS, A1 = Allele1, A2 = Allele2, original_BETA = BETA, original_SE = SE, original_VAR = var, original_PVAL = p.value) 

all_snps = cases_vs_control_pval %>% filter(sig == TRUE) %>% select(Chromosome = CHR, Position = BP, A1 = REF, A2 = ALT, rsid = ID, 
                                                                    CureGN_LRT_PVAL = P) %>% mutate(Chromosome = str_remove(Chromosome, "chr"))
allsnps_merge = merge(all_snps, sumstat_filter, 
                      by = c("Chromosome", "Position", "A1", "A2")) %>% 
  select(Chromosome, Position, A1, A2, rsid, CureGN_LRT_PVAL, original_PVAL, everything())

allsnps_merge_anno = merge(allsnps_merge, merge_nexus, by = c("Chromosome", "Position", "A1", "A2", "rsid")) %>% 
  select(Chromosome, Position, A1, A2, rsid, everything()) %>% arrange(original_PVAL)
write.table(allsnps_merge_anno, paste0(dir, "curegn/AF/annotated_results/curegn_replication_ns_allsig.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")
