library(GMMAT)
library(genio)
library(tidyverse)
library(SNPRelate)

# prep SMMAT input
# extract group file for european data
bim <- read_bim("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20.bim")
bim_clean = bim %>% select(Chr = chr, Otherinfo5 = pos, Alt = alt, Ref = ref, id)
gene_set = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/SKAT/gene_group_0321.txt", sep = " ", header = TRUE)
group_file = merge(bim_clean, gene_set, by = c("Chr", "Otherinfo5", "Ref", "Alt"), all.x = TRUE) %>% drop_na(Gene.wgEncodeGencodeBasicV26) %>% mutate(weight = 1) %>% select(Gene.wgEncodeGencodeBasicV26, Chr, Otherinfo5, Ref, Alt, weight)
group_file = group_file[order(group_file$Gene.wgEncodeGencodeBasicV26),]
print(length(unique(group_file$Gene.wgEncodeGencodeBasicV26)))
print(is.unsorted(group_file$Gene.wgEncodeGencodeBasicV26))
group.file <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr//SKAT/group_file_allPC20.txt"
write.table(group_file, group.file, row.names=FALSE, quote=FALSE, col.names=FALSE,sep = "\t")
# variant_list = group_file %>% select(id)
# write.table(variant_list, "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/SKAT/smmat_w_subset.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

# prep model equation:
name <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20.grm"
grm = as.matrix(read_grm(name)$kinship)

phenofile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_pheno.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) %>% mutate(V3 = ifelse(V3 == 1, 0, 1))
covarfile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_covar.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 
eigenfile = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20.eigenvec", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 
pheno = merge(phenofile %>% dplyr::rename(trait = V3), covarfile %>% dplyr::rename(sex = V3, race = V4), by = c("V1","V2")) 
pheno_all = merge(pheno, eigenfile, by = c("V1","V2")) %>% dplyr::rename(famid = V1, id = V2)
model0 <- glmmkin(trait ~  sex + race + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, data = pheno_all, kins = grm, id = "id",family = binomial(link = "logit"))


geno.file <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/SKAT/allssns_ctr_mac20_seq.gds"

print("start SMMAT")
# SMMAT
out <- SMMAT(model0, group.file = group.file,  geno.file = geno.file, MAF.range = c(1e-7, 0.05), miss.cutoff = 1, method = "davies", tests = "S", use.minor.allele = TRUE)
write.table(out, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/SKAT/SMMAT_out/allssns_PC20_smmat_sameweight_maxaf005_0424.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)
print("SMMAT done")
