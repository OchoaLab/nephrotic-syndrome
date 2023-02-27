### extract bristol data
# plink2 --bfile /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas.all.b38.n2656/nephrotic.syndrome.gwas.all.b38.n2656.R2 --keep /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/bristol_remove.txt --out NS_Bristol --make-bed

### subset SNPs from bigger data (NS + TGP) - need to recode
library(tidyverse)
library(genio)
bim <- read_bim("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/ssns_tpg_merge.bim")
# bim_id = bim %>%
#   mutate(SNP = paste0(chr, ":", pos, "-", ref, "-", alt)) %>% select(SNP)
write.table(bim %>% select(id), "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/NS_TGP_SNPlist.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)


bim_b = read_bim("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/NS_Bristol.bim")
bim_recode = bim_b %>% mutate(id = paste0(chr, ":", pos))

bristol = read_plink("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/NS_Bristol.bed")
bristol_X = bristol$X
bristol_X[1:5, 1:5]
rownames(bristol_X) <- bim_recode$id
write_plink( "/datacommons/ochoalab/ssns_gwas/replication/bristol_data/NS_Bristol_recode", bristol_X, bim_recode, bristol$fam )

## dedup after subset
bim_sub = read_bim("/datacommons/ochoalab/ssns_gwas/replication/bristol_data/NS_Bristol_subset.bim")

bim_id = bim_sub %>% group_by(id) %>% filter(n()>1) %>% pull(id) 

ind_remove = c()
for (id in unique(bim_id)) {
  
  ind = grep(id, rownames(X))
  x_temp = X[ind,]
  print(mean(x_temp[1,] == x_temp[2,], na.rm = TRUE))
  if (mean(x_temp[1,] == x_temp[2,], na.rm = TRUE) > 0.99) {
    ind_remove = c(ind_remove, ind[2])
  } else {
    ind_remove = c(ind_remove, ind)
  }
}

X <- X[ -ind_remove, ]
bim <- bim[ -ind_remove, ]
write_plink( "NS_Bristol_dedup", X, bim, data$fam )

### flip necessary SNPs
## flip 
#time plink --bfile NS_Bristol_dedup --flip /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/flip_10.txt --snps-only just-acgt --make-bed --out NS_Bristol_flip
## remove snps
#time plink --bfile NS_Bristol_flip --exclude /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/remove_10_new.txt --make-bed --out NS_Bristol_flip_remove

### run genotype missingness filer
## time plink2 --bfile NS_Bristol_flip_remove --mind 0.1 --make-bed --out NS_Bristol_mind

### for imputation: create vcf file and split by chromosome
# module load Plink/1.90
# time plink --bfile NS_Bristol_flip_remove --recode vcf --out NS_Bristol
# module unload Plink/1.90
# module load R/3.6.0
# module load bcftools/1.4
# module load Plink/2.00a2LM
# for i in {1..22}
# do
# 
# plink2 --bfile ../NS_Bristol_flip_remove --chr ${i} --output-chr chrM --export vcf bgz id-paste=iid --out chr_${i}
# bcftools index chr_${i}.vcf.gz
# 
# done
# 
# module unload R/3.6.0
# module unload bcftools/1.4
# module unload Plink/2.00a2LM

### after imputation:
###concatenate all the separated vcf files by chromosome into one file
# module load bcftools/1.4
# bcftools concat -O z -o bristol_impute.vcf.gz chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz
# module unload bcftools/1.4