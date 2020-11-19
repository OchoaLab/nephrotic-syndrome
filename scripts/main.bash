# steps for intersecting loci between GWAS and TGP datasets, for merging

# go where the GWAS data is
cd ~/dbs/ssns_gwas/

wc -l ssns.gwas.hg37.{bim,fam}
# 1748250 ssns.gwas.hg37.bim
#     762 ssns.gwas.hg37.fam

perl -w common_loci_bim.pl ssns.gwas.hg37.bim ~/dbs/tgp/plink2/all_phase3_filt-minimal.bim out.bim
# Reading file 1: ssns.gwas.hg37.bim
# File 1 (ssns.gwas.hg37.bim): 1677135 loci passed filters
# Reading file 2: /home/viiia/dbs/tgp/plink2/all_phase3_filt-minimal.bim
# and writing: out.bim

wc -l out.bim 
# 1407537 out.bim
gzip out.bim

################################

# script for plotting allele frequencies in cases vs controls
# here we provide name of BED data in DMPI server, but analysis will work with other BED files as long as case/control status is given in FAM pheno column as 2/1 values.
Rscript af-case-ctrl.R /dmpi/analysis/Rasheed/ssns_gwas/merge_tgp/combine
