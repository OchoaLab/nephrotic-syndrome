library(tidyverse)
library(qqman)
output = read.table('/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/saige_output.txt', header = TRUE)
#head(output)
#output %>% filter(POS == 199818960)
output_clean = output %>% select(Chromosome = CHR, Position = POS, A1 = Allele1, A2 = Allele2, BETA, SE, Tstat, VAR = var, PVAL = p.value,
                                 AF_case, N_case, AF_ctrl, N_ctrl) %>% filter(PVAL < 1e-5)

hist(output$p.value, main = "srns vs ctrl (afr)")
hist(output$p.value.NA)
qq(output$p.value, main = "srns vs ctrl (afr)")
man_plot = output %>% select(SNP = MarkerID, CHR, BP = POS, P = p.value)
manhattan(man_plot)

# upload to locus zoom
results = output %>% select(SNP = MarkerID, CHR, POS, P = p.value, A1 = Allele1, A2 = Allele2) %>% drop_na(P)
write.table(results, #%>% arrange(P) %>% head(10000) %>% arrange(CHR, BP),
            "/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/afr/LZ_srns_ctrl_afr.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")


# upload to snpnexus
SNP_nexus = output %>% filter(p.value < 1e-5) %>% select(Id = CHR, Position = POS, Allele1, Allele2) %>% 
  #mutate(CHR = str_remove(CHR, "chr")) %>% 
  mutate(type = "Chromosome", Strand = 1) %>% 
  select(type, Id, Position, Allele1, Allele2, Strand) %>% distinct()
write.table(SNP_nexus,  "/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/afr/snpnexus_srns_ctrl_afr.txt" , 
            row.names=FALSE, quote=FALSE, col.names = FALSE, sep = "\t")

neargens = read.table("/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/afr/near_gens.txt", header = TRUE, sep = "\t") %>% select(-Chromosome, -Position)
gencoord = read.table("/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/afr/gen_coords.txt", header = TRUE, sep = "\t") %>% 
  select(Variation.ID, Chromosome, Position, A1 = REF.Allele, A2 = ALT.Allele..IUPAC., rsid = dbSNP)
merge_nexus = merge(gencoord, neargens, by = "Variation.ID") %>% select(-Variation.ID) %>% 
  mutate(A1 = ifelse(A1 == TRUE, "T", A1))


# annotated clean data
annotated_df = merge(output_clean, merge_nexus, by = c("Chromosome", "Position", "A1", "A2")) %>% 
  select(Chromosome, Position, A1, A2, rsid, everything()) %>% arrange(PVAL)
write.csv(annotated_df, "/datacommons/ochoalab/ssns_gwas/saige/annotated/srns_ctrl_afr.csv", quote = FALSE)
