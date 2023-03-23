### annotate data using ANNOVAR first, requires vcf file

## converts vcf to "avinput" file, takes a long time
# perl table_annovar.pl ../all_ssns_ctr_hwe_mac20.vcf tempdir/ -buildver hg38 -out all_ssns_anno -remove -protocol wgEncodeGencodeBasicV26 -operation g -nastring . -vcfinput -polish --verbose
# awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}' all_ssns_anno.avinput > all_ssns_anno_cols.avinput
# perl table_annovar.pl all_ssns_anno_cols.avinput dir/ -buildver hg38 -outfile all_ssns_anno_step2_0215 -remove -protocol wgEncodeGencodeBasicV26 -operation g -nastring . -polish --verbose -otherinfo

### extract gene group info we need for SMMAT/SKAT analysis
file = read.table("all_ssns_anno_step2_0215.hg38_multianno.txt", sep = "\t", header = TRUE)
file_filter = file %>% 
  filter(Func.wgEncodeGencodeBasicV26 == "exonic" | Func.wgEncodeGencodeBasicV26 == "exonic;splicing") %>% 
  filter(ExonicFunc.wgEncodeGencodeBasicV26 != "synonymous SNV") %>% 
  select(Chr, Otherinfo5, Ref, Alt, Func.wgEncodeGencodeBasicV26, Gene.wgEncodeGencodeBasicV26)

gene_group = separate_rows(file_filter, Gene.wgEncodeGencodeBasicV26 ,sep = ";") %>% select(-Func.wgEncodeGencodeBasicV26) %>% distinct(.keep_all = TRUE)

length(unique(gene_group$Gene.wgEncodeGencodeBasicV26))

write.table(gene_group, "gene_group.txt", row.names=FALSE, quote=FALSE, col.names=TRUE)

### SMMAT requires gds file
bed.fn <- "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr_hwe_mac20.bed"
fam.fn <- "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr_hwe_mac20.fam"
bim.fn <- "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr_hwe_mac20.bim"
library(SeqArray)
seqBED2GDS(bed.fn, fam.fn, bim.fn, "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/SKAT/allssns_mac20_seq.gds")

### run skat_allPC_20.R file, ~ 12 hours


### read output file from SMMAT
all_skat <- read.table("allssns_PC20_smmat_sameweight.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE)
names(all_skat) = c("variant", "variant_N", "min_missing", "mean_missing", "max_missing", "min_AF", "mean_AF", "max_AF", "SKAT_PVAL")
write_df = all_skat %>% arrange(SKAT_PVAL) %>% filter(SKAT_PVAL < 3e-6)
write.table(write_df, "variantanalysis_0220.txt", row.names=FALSE, quote=FALSE, col.names=TRUE, sep = "\t")