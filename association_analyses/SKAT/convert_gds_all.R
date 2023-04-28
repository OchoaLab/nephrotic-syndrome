library(SNPRelate)
bed.fn <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20.bed"
fam.fn <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20.fam"
bim.fn <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20.bim"

# convert
#snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/SKAT/w_mac20.gds")

library(SeqArray)
seqBED2GDS(bed.fn, fam.fn, bim.fn, "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/SKAT/allssns_ctr_mac20_seq.gds")
