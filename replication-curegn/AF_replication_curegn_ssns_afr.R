library(tidyverse)
library(genio)
library(qqman)
library(popgeninfer)
dir = '/datacommons/ochoalab/ssns_gwas/replication/saige_results/'
setwd( dir )

#############################
disease_subtype = "ssns_ctrl"
ancestry = "afr"
#############################

# step 1: extract gnomad SNP id's
gnomad_AC <- read.table(paste0(dir, "gnomad/outfile/gnomad-genome_ssns_ctrl_afr_clean.txt"), sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
gnomad_snp = gnomad_AC %>% mutate(SNP = paste0(CHROM, ":", POS)) %>% mutate(SNP = str_remove(SNP, "chr")) %>% select(SNP)
write.table(gnomad_snp, paste0(dir, "curegn/AF/ssns/afr/gnomad_snp_id.txt"),
            row.names=FALSE, quote=FALSE, col.names = FALSE)

# test = read_bim("/datacommons/ochoalab/curegn/merge_tgp/admixture/curegn_NS_afr.bim")
# curegn only has cases, no controls, so skip step 2 of removing imputation biased snps from control vs control
# run LRT on CureGN cases vs Gnomad controls
# run allele_freq_curegn.q, process outfiles

afr_ac <- read.table(paste0(dir, "curegn/AF/ssns/afr/curegnAC_afr_sub.acount")) 
colnames(afr_ac) = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
afr_ac_clean = afr_ac  %>%
  separate(ID, c("CHROM", "POS"))  %>% dplyr::rename(curegn_AC_afr = ALT_FREQS, curegn_AN_afr = OBS_CT)

cases_vs_control = merge(gnomad_AC, afr_ac_clean %>% mutate(CHROM = paste0("chr", CHROM)), by = c("CHROM", "POS", "REF", "ALT")) %>% distinct()

# LRT on curegn vs Gnomad (case vs control)
# define function input
# x1/n1: curegn; x2/n2: gnomad
x1 <- cbind(cases_vs_control$curegn_AC_afr)
n1 <- cbind(cases_vs_control$curegn_AN_afr)
x2 <- cbind(cases_vs_control$AC_afr)
n2 <- cbind(cases_vs_control$AN_afr)
# p-values for forward alignment (only alignment for non-revcomp SNPs)
pvals_b <- af_test( x1, n1, x2, n2 )$pval
# define p-value threshold
samplesize = length(pvals_b)
pcut = 0.05/samplesize
cases_vs_control_test = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% 
  dplyr::rename(CHR = CHROM, BP = POS, P = pvals_b) %>% mutate(CHR = as.numeric(str_remove(CHR, "chr")), SNP = paste0(CHR, ":", BP)) %>% 
  arrange(CHR, BP) %>% mutate(sig = ifelse(P < pcut, TRUE, FALSE))

output_plot = cbind(cases_vs_control, pvals_b)
write.table(output_plot, paste0(dir, "curegn/AF/annotated_results/replication_pvals_", disease_subtype,"_afr.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")

# define p-value threshold
unclumped_sig = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < pcut) 
nrow(unclumped_sig)

hist(pvals_b, breaks = 50, main = "curegn SSNS vs gnomad controls: AFR")
qq(pvals_b, main = "curegn SSNS vs gnomad controls: AFR")

### write file for LD clump 
curegn_snps = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% mutate(CHROM = str_remove(CHROM, "chr")) %>% 
  mutate(SNP = paste0(CHROM, ":", POS)) %>% select(SNP, P = pvals_b) %>% distinct()
write.table(curegn_snps, paste0(dir, "curegn/AF/ssns/afr/replication_unclump_curegn_snps_afr.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")
### clump p-values to get effective number of tests
clump_curegn <- read.table(paste0(dir,  "curegn/AF/ssns/afr/clump_curegn_1233.clumped"), header = TRUE) 
nrow(clump_curegn)

# apply new p-value threshold
pval_filter = cbind(cases_vs_control, pvals_b) %>% arrange(pvals_b) %>% filter(pvals_b < 0.05/nrow(clump_curegn)) %>% 
  mutate(SNP = paste0(CHROM, ":", POS, ":", REF, ":", ALT))
nrow(pval_filter)
length(intersect(clump_curegn$SNP, pval_filter$SNP))
intersect_list = intersect(clump_curegn$SNP, pval_filter$SNP)
independent_regions = pval_filter %>% filter(SNP %in% intersect_list) %>% select(CHROM, POS, REF, ALT, ID, SNP, pvals_b) %>% 
  mutate(CHROM = str_remove(CHROM, "chr"))


# include summary stats in output 
ancestry_all = read.table(paste0("/datacommons/ochoalab/ssns_gwas/saige/ssns_ctrl/afr/saige_output.txt" ),  sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
# annotate snps according to ld clumped number of independent snps

cases_vs_control_pval = cases_vs_control_test %>% mutate(sig = ifelse(P < (0.05/nrow(clump_curegn)), TRUE, FALSE))
table(cases_vs_control_pval$sig)

ssns_repli_table = cases_vs_control_test %>% filter(ID == "rs57588792")
gnomad_AC %>% filter(ID == "rs57588792")
afr_ac  %>% filter(ID == "rs57588792")

SNP_nexus = cases_vs_control_pval %>% filter(sig == TRUE) %>% select(CHR, BP, REF, ALT, ID, SNP, P) %>% 
  mutate(type = "Chromosome", Strand = 1) %>% 
  select(type, Id = CHR, Position = BP, Allele1 = REF, Allele2 = ALT, Strand) %>% distinct()
write.table(SNP_nexus, paste0(dir,"curegn/AF/ssns/afr/curegn_snp_nexus_allsig_afr.txt"), 
            row.names=FALSE, quote=FALSE, col.names = FALSE, sep = "\t")

neargens = read.table(paste0(dir,"curegn/AF/ssns/afr/near_gens.txt"), header = TRUE, sep = "\t") %>% select(-Chromosome, -Position)
gencoord = read.table(paste0(dir, "curegn/AF/ssns/afr/gen_coords.txt"), header = TRUE, sep = "\t") %>% 
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

allsnps_merge_anno = merge(allsnps_merge, merge_nexus, by = c("Chromosome", "Position", "A1", "A2", "rsid"), all.x = TRUE) %>% 
  select(Chromosome, Position, A1, A2, rsid, everything()) %>% arrange(original_PVAL) %>% distinct()

write.table(allsnps_merge_anno, paste0(dir, "curegn/AF/annotated_results/curegn_replication_ssns_allsig_afr.txt"), 
            row.names=FALSE, quote=FALSE, col.names = TRUE, sep = "\t")
