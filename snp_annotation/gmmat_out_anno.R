library(tidyverse)

# Annotate gmmat output using SNPNexus: https://www.snp-nexus.org/
# annotate with suggestive p-value threshold (P < 1 × 10−5)

setwd('/datacommons/ochoalab/ssns_gwas/imputed/')

#############################
disease_subtype = "ns_ctrl"
ancestry = "all"
#############################
if (disease_subtype == "ns_ctrl" & ancestry == "all"){
  data = read.table("mac20-glmm-score.txt", header = TRUE) 
} else {
  data = read.table(paste0(ancestry, "/mac20-glmm-score.txt"), header = TRUE) 
}

if (disease_subtype != "ns_ctrl" & ancestry == "all"){
  data = read.table(paste0(disease_subtype, "/mac20-glmm-score.txt"), header = TRUE) 
  } else {
  data = read.table(paste0(disease_subtype, "/", ancestry, "/mac20-glmm-score.txt"), header = TRUE) 
}

data_subset = data %>% filter(PVAL < 1e-05) %>% mutate(CHROM = "Chromosome", STRAND = 1) %>% select(CHROM, CHR, POS, REF, ALT, STRAND)
## file to upload to SNPNexus
write.table(data_subset, paste0("annotate/", disease_subtype, "_", ancestry, "_nexus.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
#############################
## read SNPNexus files
neargens = read.table(paste0("annotate/near_gens_", disease_subtype, "_", ancestry, ".txt"), header = TRUE, sep = "\t") %>% select(-Chromosome, -Position)
gencoord = read.table(paste0("annotate/gen_coords_", disease_subtype, "_", ancestry, ".txt"), header = TRUE, sep = "\t") %>% 
  select(Variation.ID, Chromosome, Position, A1 = REF.Allele, A2 = ALT.Allele..IUPAC., rsid = dbSNP)
merge_nexus = merge(gencoord, neargens, by = "Variation.ID") %>% select(-Variation.ID)
#A2 is the minor allele?
# subset gmmat summary stats
gmmat_data = data %>% filter(PVAL < 1e-05) %>% select(SNP, CHR, POS, A1 = REF, A2 = ALT, SCORE, VAR, PVAL) 
# get case control AF
## write subset snp id for plink
write.table(gmmat_data %>% select(SNP),
            paste0("annotate/snp_id/", disease_subtype, "_", ancestry, "_snpid.txt"), row.names=FALSE, quote=FALSE, col.names = FALSE)
##############################
## read af case control outut
af_cc = read.table(paste0("annotate/", disease_subtype, "_", ancestry, "_casecontrol.frq.cc"), header = TRUE) %>% 
  select(SNP, CHR, A1, A2, N_case_allele = NCHROBS_A, AF_case = MAF_A, N_control_allele = NCHROBS_U, AF_control = MAF_U)

merge_sumstat = merge(merge_nexus, gmmat_data, by.x = c("Chromosome", "Position", "A1", "A2"), 
                      by.y = c("CHR", "POS", "A1", "A2"), all = TRUE) %>% select(-A1, -A2)
# a1/a2 are flipped from af calculations
# plink: a1 is alternate; --freq calculates minor allele freq (alt)
merge_af = merge(merge_sumstat, af_cc, by = "SNP") %>% select(Chromosome, Position, A1, A2, everything(), -CHR, -SNP) %>% arrange(PVAL)
write.csv(merge_af, paste0("annotate/clean/", disease_subtype, "_", ancestry, ".csv"), row.names = FALSE)


# write file for locuszoom
write.table(gmmat_data, paste0("annotate/locuszoom/", disease_subtype, "_", ancestry, ".txt"), sep = "\t", row.names=FALSE, quote=FALSE, col.names = TRUE)
